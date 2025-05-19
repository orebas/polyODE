#!/usr/bin/env python3
import json, sys

# Read PHC output from stdin
text = sys.stdin.read().strip()

try:
    # Directly evaluate the Python dictionary output from PHC
    data = eval(
        text, {"__builtins__": {}}, {}
    )  # Disable built-in functions for security
except Exception as e:
    print(f"Eval Failed: {e}", file=sys.stderr)
    print("Problematic Input:", file=sys.stderr)
    print(text[:1000], file=sys.stderr)
    sys.exit(1)


# Convert complex numbers into JSON-friendly {"re": x, "im": y}
def fix_complex(obj):
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    elif isinstance(obj, dict):
        return {k: fix_complex(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [fix_complex(x) for x in obj]
    return obj


data_fixed = fix_complex(data)
json.dump(data_fixed, sys.stdout, indent=2) 