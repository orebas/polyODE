#!/usr/bin/env python3
import sys
import json
import re
import os

def parse_rational_string(s_val_raw):
    s_val = str(s_val_raw).strip()
    if '/' in s_val:
        parts = s_val.split('/')
        num_str = parts[0].strip()
        den_str = parts[1].strip()
        try:
            num = float(num_str)
            if '^' in den_str:
                den_parts = den_str.split('^')
                base_str = den_parts[0].strip()
                exp_str = den_parts[1].strip()
                if not base_str or not exp_str:
                    return None
                base = float(base_str)
                exponent = float(exp_str)
                den = pow(base, exponent)
            else:
                if not den_str:
                    return None
                den = float(den_str)
            if den == 0:
                return None
            return num / den
        except ValueError:
            return None
        except Exception: # Catch any other error during parsing
            return None
    else:
        try:
            return float(s_val)
        except ValueError:
            return None

def extract_solution_points_list_str(content_cleaned):
    # Updated regex: r'\[\s*0\s*,\s*\[\s*1\s*,\s*(\[.*\])\s*\]\s*\]' to capture the entire solutions block.
    # The (.*) is greedy and DOTALL allows . to match newlines.
    match_main_structure = re.search(r'\[\s*0\s*,\s*\[\s*1\s*,\s*(\[.*\])\s*\]\s*\]', content_cleaned, re.DOTALL)
    if match_main_structure:
        solutions_block = match_main_structure.group(1).strip()
        return solutions_block
    return None

def parse_msolve_output(file_path):
    try:
        with open(file_path, 'r') as f:
            content = f.read()
    except Exception as e:
        return {"status": "error", "message": f"Failed to read msolve output file: {str(e)}"}

    content_cleaned = content.strip()
    if content_cleaned.endswith(':'):
        content_cleaned = content_cleaned[:-1]

    if content_cleaned == "[-1]":
        return {"status": "no_solution", "message": "Inconsistent system."}
    if re.fullmatch(r'\[1\s*,\s*\d+\s*,\s*-1\s*,\s*\[\]\s*\]', content_cleaned): # Check for positive dimensional / infinite solutions
        return {"status": "infinite_solutions", "message": "Positive dimensional system."}

    solutions_list_block_str = extract_solution_points_list_str(content_cleaned)

    if solutions_list_block_str is None:
        return {"status": "unknown_format", "message": f"Could not extract solutions block from: {content_cleaned[:200]}"}
    if solutions_list_block_str == "[]":
        return {"status": "no_real_roots", "message": "No real solutions found (solutions list is '[]')."}

    solutions_data = []
    
    # Primary strategy: Use regex to find all individual [low, high] pairs
    # This regex finds the content within the innermost brackets representing value intervals.
    # It's made non-greedy for the content capture to handle adjacent pairs correctly.
    value_interval_strs = re.findall(r'\[\s*([^\[\]]+?)\s*\]', solutions_list_block_str)

    parsed_values_from_intervals = []
    for interval_content_str in value_interval_strs:
        parts = interval_content_str.split(',')
        if len(parts) == 2:
            low_str = parts[0].strip()
            high_str = parts[1].strip()
            
            low_val = parse_rational_string(low_str)
            high_val = parse_rational_string(high_str)

            if low_val is not None and high_val is not None:
                parsed_values_from_intervals.append((low_val + high_val) / 2.0)
            else:
                # If any rational string parsing fails, this path is problematic for the primary strategy
                parsed_values_from_intervals = [] # Clear to indicate failure of this strategy
                break 
        else:
            # Malformed pair content
            parsed_values_from_intervals = [] # Clear
            break
    
    if parsed_values_from_intervals:
        # Heuristic to differentiate between:
        # 1. A single solution point for multiple variables (e.g., "[[[v1_l,v1_h],[v2_l,v2_h]]]")
        # 2. Multiple solution points, each for a single variable (e.g., "[[[v1_l,v1_h]],[[v2_l,v2_h]]]")
        
        # If the entire block is a list containing lists of pairs, it's likely case 1.
        # A simple check: count the depth of brackets at the start.
        # More robustly: if the number of value_interval_strs matches the number of items after json.loads of the outer structure.
        
        # If the solutions_list_block_str starts with "[[[" and ends with "]]]" and does NOT contain the pattern
        # "]], [[" (which separates multiple solution points), then it's likely a single solution point with
        # multiple variables. Otherwise, treat each parsed value as its own solution point (single-variable system
        # with multiple distinct roots).
        if solutions_list_block_str.startswith('[[[') and solutions_list_block_str.endswith(']]]') and \
           re.search(r'\]\],\s*\[\[', solutions_list_block_str) is None:
            solutions_data.append(parsed_values_from_intervals)
        else:
            # Otherwise, assume each parsed value is a separate solution point (for a single variable system)
            for val in parsed_values_from_intervals:
                solutions_data.append([val])
    
    # Fallback strategy: If regex strategy failed or yielded no data, try json.loads on the cleaned block
    # This is particularly for cases where msolve might output pure numbers in a way that's valid JSON
    # but not caught by the interval regex (e.g. if msolve output single numbers not in [low,high] pairs)
    if not solutions_data and solutions_list_block_str:
        solutions_list_block_str_cleaned_for_json = solutions_list_block_str.replace('\\n', ' ')
        try:
            # print(f"DEBUG_PY: Attempting fallback json.loads on: '{solutions_list_block_str_cleaned_for_json}'", file=sys.stderr)
            msolve_parsed_solutions_list = json.loads(solutions_list_block_str_cleaned_for_json)
            if isinstance(msolve_parsed_solutions_list, list):
                # Assume the structure is already what C++ expects: list of lists of numbers
                # This part needs to be careful if msolve_parsed_solutions_list is like [[0.5, 0.75]] for a 2-var system
                # or like [[0.5], [0.75]] for a 1-var system with two roots.
                # The C++ side expects each inner list to be a full solution point.
                
                # Check if it's a list of numbers (single solution point for multi-var)
                # or a list of lists of numbers (multiple solution points)
                if all(isinstance(item, (int, float)) for item in msolve_parsed_solutions_list):
                    # It's a flat list like [0.5, 0.75], treat as single solution
                    solutions_data.append(msolve_parsed_solutions_list)
                elif all(isinstance(sublist, list) for sublist in msolve_parsed_solutions_list):
                    # It's already a list of lists like [[0.5], [0.75]] or [[0.5, 0.75]]
                    solutions_data = msolve_parsed_solutions_list
        except json.JSONDecodeError:
            # print(f"DEBUG_PY: Fallback json.loads failed for: '{solutions_list_block_str_cleaned_for_json}'", file=sys.stderr)
            pass # Fallback failed, do nothing, rely on primary strategy's outcome or lack thereof

    if solutions_data:
        return {"status": "success", "solutions": solutions_data}
    else: # If both primary and fallback strategies yield nothing
        return {"status": "unknown_format", "message": f"Could not parse solutions from block: {solutions_list_block_str[:200]}"}

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)
    file_path = sys.argv[1]
    result = parse_msolve_output(file_path)
    print(json.dumps(result))