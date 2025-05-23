# Contributing to PolyODE

Thank you for your interest in contributing to PolyODE! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/yourusername/polyODE.git`
3. Follow the build instructions in [BUILDING.md](BUILDING.md)
4. Create a feature branch: `git checkout -b feature/your-feature-name`

## Development Setup

### Prerequisites
- Follow the setup instructions in [BUILDING.md](BUILDING.md)
- Ensure all tests pass: `cd build && ctest`

### Code Style
- Use C++17 features where appropriate
- Follow existing naming conventions
- Include header guards and proper includes
- Document public APIs with comments

### Testing
- Add tests for new functionality in the `tests/` directory
- Ensure all existing tests pass
- Use Google Test framework for unit tests
- Test files should be named `*_test.cpp`

## Contribution Process

### 1. Issue Discussion
- Check existing issues before creating new ones
- Discuss significant changes in issues before implementation
- Tag issues appropriately (bug, feature, documentation, etc.)

### 2. Code Changes
- Make focused, logical commits
- Write clear commit messages
- Keep changes as small as practical
- Update documentation as needed

### 3. Pull Request
- Create a pull request against the `main` branch
- Include a clear description of changes
- Reference any related issues
- Ensure CI checks pass

## Areas for Contribution

### High Priority
- Bug fixes and stability improvements
- Performance optimizations
- Documentation improvements
- Additional test coverage

### Features
- New ODE system examples
- Parameter estimation algorithms
- Identifiability analysis methods
- Visualization tools

### Infrastructure
- Build system improvements
- CI/CD enhancements
- Cross-platform compatibility
- Package management

## Code Guidelines

### Header Files
```cpp
#pragma once

#include <standard_library>
#include "local_header.hpp"

namespace polyode {
    // Your code here
}
```

### Implementation Files
```cpp
#include "header.hpp"

#include <standard_library>
#include "other_local_headers.hpp"

namespace polyode {
    // Implementation
}
```

### Testing
```cpp
#include <gtest/gtest.h>
#include "your_header.hpp"

TEST(TestSuiteName, TestCaseName) {
    // Test implementation
    EXPECT_EQ(expected, actual);
}
```

## Documentation

### Code Documentation
- Document public APIs with clear comments
- Include parameter descriptions and return values
- Provide usage examples for complex functions

### README Updates
- Update README.md for new features
- Add examples to demonstrate functionality
- Update dependency lists if needed

## Getting Help

- Open an issue for questions
- Check existing documentation first
- Provide minimal reproducible examples for bugs

## License

By contributing, you agree that your contributions will be licensed under the same license as the project.

## Recognition

Contributors will be acknowledged in the project documentation and releases.

Thank you for contributing to PolyODE!