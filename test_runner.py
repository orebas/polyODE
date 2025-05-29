#!/usr/bin/env python3
"""
Run Google Tests individually with timeout and error handling in a CMake/CTest project.
Can work with either CTest or direct test binaries.
"""

import subprocess
import sys
import time
import signal
import argparse
import json
import re
import os
from pathlib import Path
from collections import defaultdict

class TestResult:
    PASSED = "PASSED"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    TIMEOUT = "TIMEOUT"
    SEGFAULT = "SEGFAULT"
    ERROR = "ERROR"

class CTestRunner:
    def __init__(self, build_dir=".", timeout=60, verbose=False):
        self.build_dir = Path(build_dir)
        self.timeout = timeout
        self.verbose = verbose
        self.results = defaultdict(list)
        
    def find_test_binaries(self):
        """Find all test binaries in the build directory"""
        test_binaries = []
        
        # Look for files ending with _test or _tests
        for path in self.build_dir.iterdir():
            if path.is_file() and os.access(path, os.X_OK):
                # Check if filename ends with _test or _tests
                if path.name.endswith('_test') or path.name.endswith('_tests'):
                    # Verify it's a Google Test binary
                    try:
                        result = subprocess.run(
                            [str(path), '--gtest_list_tests'],
                            capture_output=True,
                            timeout=5
                        )
                        if result.returncode == 0:
                            test_binaries.append(path)
                            if self.verbose:
                                print(f"Found test binary: {path.name}")
                    except Exception as e:
                        if self.verbose:
                            print(f"Skipping {path.name}: {e}")
                        
        # Also search subdirectories if needed
        for subdir in self.build_dir.iterdir():
            if subdir.is_dir() and subdir.name not in ['CMakeFiles', 'Testing', 'vcpkg_installed']:
                for path in subdir.glob('*_test*'):
                    if path.is_file() and os.access(path, os.X_OK):
                        try:
                            result = subprocess.run(
                                [str(path), '--gtest_list_tests'],
                                capture_output=True,
                                timeout=5
                            )
                            if result.returncode == 0:
                                test_binaries.append(path)
                        except:
                            pass
                        
        return sorted(test_binaries)
    
    def get_ctest_tests(self):
        """Get all tests registered with CTest"""
        try:
            # Run ctest with -N to list tests without running them
            result = subprocess.run(
                ['ctest', '-N'],
                cwd=self.build_dir,
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                print(f"Error listing CTest tests: {result.stderr}")
                return []
            
            # Parse the output to get test names
            tests = []
            for line in result.stdout.splitlines():
                # Look for lines like "  Test #1: test_name"
                match = re.match(r'\s*Test\s+#\d+:\s+(.+)', line)
                if match:
                    tests.append(match.group(1))
                    
            return tests
            
        except FileNotFoundError:
            print("ctest not found. Make sure CMake is installed.")
            return []
        except Exception as e:
            print(f"Error listing tests: {e}")
            return []
    
    def run_ctest_single(self, test_name):
        """Run a single test through CTest"""
        if self.verbose:
            print(f"Running (via ctest): {test_name}")
            
        start_time = time.time()
        
        try:
            # Run the test with CTest
            result = subprocess.run(
                ['ctest', '-R', f'^{re.escape(test_name)}$', '-V'],
                cwd=self.build_dir,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            duration = time.time() - start_time
            
            # Parse CTest output
            if "Passed" in result.stdout and result.returncode == 0:
                return TestResult.PASSED, duration, ""
            elif "***Failed" in result.stdout or result.returncode != 0:
                # Check for segfault indicators
                if "Segmentation fault" in result.stdout or "signal 11" in result.stdout:
                    return TestResult.SEGFAULT, duration, "Segmentation fault"
                else:
                    return TestResult.FAILED, duration, result.stdout
            else:
                return TestResult.ERROR, duration, "Unknown test result"
                
        except subprocess.TimeoutExpired:
            # Kill the test process
            subprocess.run(['pkill', '-f', test_name], capture_output=True)
            return TestResult.TIMEOUT, self.timeout, f"Test timed out after {self.timeout} seconds"
            
        except Exception as e:
            return TestResult.ERROR, 0, f"Error running test: {str(e)}"
    
    def run_binary_tests(self, binary_path):
        """Run all tests from a single Google Test binary"""
        runner = GoogleTestRunner(str(binary_path), self.timeout, self.verbose)
        tests = runner.get_all_tests()
        
        print(f"\nRunning tests from: {binary_path}")
        print(f"Found {len(tests)} tests")
        
        for test in tests:
            status, duration, error_msg = runner.run_single_test(test)
            
            self.results[status].append({
                'name': f"{binary_path.name}::{test}",
                'duration': duration,
                'error': error_msg,
                'binary': str(binary_path)
            })
            
            if self.verbose or status not in [TestResult.PASSED, TestResult.SKIPPED]:
                self._print_test_result(f"{binary_path.name}::{test}", status, duration)
    
    def _print_test_result(self, test_name, status, duration):
        """Print a single test result"""
        status_symbol = {
            TestResult.PASSED: "✓",
            TestResult.FAILED: "✗",
            TestResult.SKIPPED: "○",
            TestResult.TIMEOUT: "⏱",
            TestResult.SEGFAULT: "💥",
            TestResult.ERROR: "⚠"
        }.get(status, "?")
        
        print(f"{status_symbol} {test_name} - {status} ({duration:.2f}s)")
    
    def run_all_tests(self, use_ctest=True):
        """Run all tests either through CTest or by finding binaries"""
        if use_ctest:
            # Try to use CTest first
            tests = self.get_ctest_tests()
            
            if tests:
                print(f"Found {len(tests)} CTest tests")
                print("-" * 60)
                
                for i, test in enumerate(tests, 1):
                    status, duration, error_msg = self.run_ctest_single(test)
                    
                    self.results[status].append({
                        'name': test,
                        'duration': duration,
                        'error': error_msg
                    })
                    
                    if self.verbose or status not in [TestResult.PASSED]:
                        print(f"[{i}/{len(tests)}] ", end="")
                        self._print_test_result(test, status, duration)
            else:
                print("No CTest tests found, searching for test binaries...")
                use_ctest = False
        
        if not use_ctest:
            # Find and run test binaries directly
            binaries = self.find_test_binaries()
            
            if not binaries:
                print("No test binaries found!")
                return
                
            print(f"Found {len(binaries)} test binaries")
            
            for binary in binaries:
                self.run_binary_tests(binary)
        
        print("-" * 60)
        self.print_summary()
    
    def print_summary(self):
        """Print summary of results"""
        total = sum(len(tests) for tests in self.results.values())
        
        print("\nSUMMARY:")
        print(f"Total tests: {total}")
        
        for status in [TestResult.PASSED, TestResult.FAILED, TestResult.SKIPPED, 
                      TestResult.TIMEOUT, TestResult.SEGFAULT, TestResult.ERROR]:
            count = len(self.results[status])
            if count > 0:
                print(f"  {status}: {count}")
                
        # Print details for problematic tests
        for status in [TestResult.FAILED, TestResult.TIMEOUT, TestResult.SEGFAULT, TestResult.ERROR]:
            if self.results[status]:
                print(f"\n{status} tests:")
                for test in self.results[status]:
                    print(f"  - {test['name']}")
                    if test.get('binary'):
                        print(f"    Binary: {test['binary']}")
                    if test['error'] and self.verbose:
                        print(f"    Error: {test['error'][:200]}...")
    
    def save_results(self, filename):
        """Save results to a JSON file"""
        with open(filename, 'w') as f:
            json.dump(dict(self.results), f, indent=2)
        print(f"\nResults saved to {filename}")

class GoogleTestRunner:
    """Runner for individual Google Test binaries"""
    def __init__(self, test_binary, timeout=20, verbose=False):
        self.test_binary = test_binary
        self.timeout = timeout
        self.verbose = verbose
        
    def get_all_tests(self):
        """Get list of all tests from the binary"""
        try:
            result = subprocess.run(
                [self.test_binary, '--gtest_list_tests'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                return []
            
            tests = []
            current_suite = None
            
            for line in result.stdout.splitlines():
                line = line.strip()
                if not line or line.startswith("Running"):
                    continue
                    
                if line.endswith('.'):
                    current_suite = line[:-1]
                elif current_suite and line:
                    test_name = line.split()[0]
                    tests.append(f"{current_suite}.{test_name}")
                    
            return tests
            
        except:
            return []
    
    def run_single_test(self, test_name):
        """Run a single test"""
        start_time = time.time()
        
        try:
            result = subprocess.run(
                [self.test_binary, f'--gtest_filter={test_name}'],
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            duration = time.time() - start_time
            
            if result.returncode == 0:
                if "[ SKIPPED ]" in result.stdout or "0 tests from" in result.stdout:
                    return TestResult.SKIPPED, duration, "Test was skipped"
                else:
                    return TestResult.PASSED, duration, ""
            elif result.returncode == -signal.SIGSEGV or result.returncode == 139:
                return TestResult.SEGFAULT, duration, "Segmentation fault"
            else:
                return TestResult.FAILED, duration, result.stderr or "Test failed"
                
        except subprocess.TimeoutExpired:
            return TestResult.TIMEOUT, self.timeout, f"Timed out after {self.timeout}s"
        except Exception as e:
            return TestResult.ERROR, 0, str(e)

def main():
    parser = argparse.ArgumentParser(description='Run tests individually with timeout handling')
    parser.add_argument('path', nargs='?', default='.', 
                       help='Build directory (for CTest) or test binary path')
    parser.add_argument('--timeout', type=int, default=20, 
                       help='Timeout per test in seconds (default: 20)')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Verbose output')
    parser.add_argument('--output', '-o', help='Save results to JSON file')
    parser.add_argument('--no-ctest', action='store_true', 
                       help='Skip CTest and search for test binaries directly')
    
    args = parser.parse_args()
    
    path = Path(args.path)
    
    if path.is_file():
        # Single test binary mode
        if not os.access(path, os.X_OK):
            print(f"Error: '{path}' is not executable")
            sys.exit(1)
            
        runner = GoogleTestRunner(str(path), args.timeout, args.verbose)
        tests = runner.get_all_tests()
        
        if not tests:
            print("No Google Tests found in binary")
            sys.exit(1)
            
        results = defaultdict(list)
        print(f"Running {len(tests)} tests from {path}")
        
        for test in tests:
            status, duration, error = runner.run_single_test(test)
            results[status].append({
                'name': test,
                'duration': duration,
                'error': error
            })
            
            if args.verbose or status not in [TestResult.PASSED]:
                print(f"{test} - {status} ({duration:.2f}s)")
                
        # Print summary
        print("\nSummary:")
        for status, tests in results.items():
            if tests:
                print(f"  {status}: {len(tests)}")
                
    else:
        # CTest/build directory mode
        runner = CTestRunner(path, args.timeout, args.verbose)
        runner.run_all_tests(use_ctest=not args.no_ctest)
        
        if args.output:
            runner.save_results(args.output)

if __name__ == "__main__":
    main()
