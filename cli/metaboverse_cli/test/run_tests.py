#!/usr/bin/env python3
"""
Run unit tests for Metaboverse redox pair mapping
"""
import unittest
import os
import sys

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the test module
from unit_tests import TestRedoxPairMapping

if __name__ == '__main__':
    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRedoxPairMapping)
    
    # Run the tests
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    
    # Print summary
    print("\n=== Test Summary ===")
    print(f"Ran {result.testsRun} tests")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    # Exit with appropriate code
    sys.exit(len(result.failures) + len(result.errors)) 