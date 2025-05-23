#!/usr/bin/env python3
import re, json, sys, os, ast
from lark import Lark, Transformer, v_args

grammar = r"""
?start: element

?element:
    | list
    | fraction
    | number
    | STRING         -> string
    | NAME

list      : "[" [element ("," element)*] "]"
fraction  : SIGNED_NUMBER "/" SIGNED_NUMBER
number    : SIGNED_NUMBER

%import common.SIGNED_NUMBER
%import common.WS
%ignore WS

// a quoted string: either "...", allowing \" escapes, or '...', allowing \' escapes
STRING  : /"(\\.|[^"\\])*"/ | /'(\\.|[^\'\\])*'/

NAME    : /[A-Za-z_]\w*/
"""

parser = Lark(grammar, parser="lalr", start="start")


@v_args(inline=True)
class TreeToJson(Transformer):
    def list(self, *items):
        return list(items)

    def fraction(self, num, den):
        return f"{num}/{den}"

    def number(self, tok):
        s = tok.value
        if re.fullmatch(r"-?\d+", s): # It's an integer string
            val = int(s)
            # Heuristic: if bit length > 62 (approximate for fitting in typical signed 64-bit after parsing)
            # or simply if it's a very long string for an integer.
            # A simpler check might be len(s) > 18 (max for signed 64-bit is 19 digits, but can be less).
            if val.bit_length() > 62: # Check bit_length for large magnitude
                return str(val) # Convert very large integers to strings for JSON
            else:
                return val # Keep smaller integers as numbers
        else: # It's a float string
            return float(s)

    def NAME(self, tok):
        return tok.value

    def string(self, tok):
        # ast.literal_eval handles both single- and double-quoted strings + escapes
        return ast.literal_eval(tok.value)


def extract_top_list(text: str) -> str:
    start = text.find("[")
    end = text.rfind("]")
    if start < 0 or end < 0:
        # If it's not a list, it might be msolve's simple error like [-1]:
        # or [1,2,-1,[]]: (positive dimensional)
        # Try to capture these simple cases as valid JSON arrays directly.
        # A more robust way would be to extend grammar, but for now, check simple list-like strings.
        cleaned_text = text.strip()
        if cleaned_text.startswith("[") and cleaned_text.endswith("]:"):
            cleaned_text = cleaned_text[:-1] # Remove trailing colon
        if cleaned_text.startswith("[") and cleaned_text.endswith("]"):
            # It might be simple enough for ast.literal_eval or json.loads directly
            # For safety, let's allow the Lark parser to try it if it's its main format
            # but if it fails the find above, let caller handle it. Or parser handles it.
            # The parser should handle valid list structures.
            # If not a list literal, it could be an msolve error output like "[-1]:"
            # The current grammar expects an element which can be a list.
            # Let's assume for now the main output we care about is a list.
            # If extract_top_list fails, the parser will fail on raw text if it's not list-like.
            # So, if it's like "[-1]:" -> becomes "[-1]", then parser makes it [-1]
            if cleaned_text.startswith("[") and cleaned_text.endswith("]"): 
                 return cleaned_text # Pass cleaned up simple list string to parser
        raise ValueError("No list literal found or unable to clean simple msolve output.")
    return text[start : end + 1]


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file> [--write-file]", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    write_to_file = False
    if len(sys.argv) > 2 and sys.argv[2] == "--write-file":
        write_to_file = True

    try:
        text = open(infile).read()
        # msolve output might have a trailing colon e.g. "[...]\n:", remove it for the parser.
        # Also handle simple cases like "[-1]:"
        cleaned_text_for_parsing = text.strip()
        if cleaned_text_for_parsing.endswith(":"):
            cleaned_text_for_parsing = cleaned_text_for_parsing[:-1].strip()
        
        # If cleaned_text_for_parsing is empty or not starting with '[', it's likely an error or unexpected format.
        if not cleaned_text_for_parsing or not cleaned_text_for_parsing.startswith("["):
            # Try to handle simple non-list outputs like msolve errors gracefully as valid JSON array of one int
            # e.g. if text was "-1" (after stripping possible "[]:")
            try:
                # Attempt to directly convert to an int and wrap in a list
                # This handles simple cases like "-1" from msolve for empty solutions
                # Needs robust check if it's really just a number
                val = int(cleaned_text_for_parsing) # check if it is an int
                data = [val] # wrap it as a list, to be JSON-like
            except ValueError: # Not an int
                # If it's not parseable as an int, and not a list, output empty JSON array
                # or a more informative error structure if preferred.
                # For now, outputting an empty list for C++ to handle as "no solution data".
                # C++ side expects a list; if it's not, nlohmann will fail.
                # So, if we can't make it a list, this will cause C++ to fail parsing, which is okay.
                # A better solution might be for this script to output {"error": "message"} json.
                # However, the Lark parser is quite flexible for list structures.
                # If it's not a list, what is it? Example: "System has dimension 1"
                # The parser currently only expects a list-like structure from extract_top_list
                # The original `extract_top_list` was safer.
                # Let's rely on the parser to fail if the structure isn't a list from extract_top_list
                pass # Let it fall through to Lark parsing attempt

        # Re-instating a safer extract_top_list call or similar logic is needed if non-list inputs are common.
        # For now, assume msolve -P 2 gives list-like output or a simple error convertible to list.
        
        # The grammar itself is defined to parse an 'element', which can be a list.
        # So, we pass the cleaned text directly to the Lark parser.
        tree = parser.parse(cleaned_text_for_parsing)
        data = TreeToJson().transform(tree)

    except Exception as e:
        print(f"Error processing file {infile}: {e}", file=sys.stderr)
        # Output an empty JSON list or an error JSON object to stdout in case of Python error
        # This allows C++ to always expect some JSON, even if it's an error indicator.
        json.dump({"python_error": str(e)}, sys.stdout, indent=2)
        sys.exit(1) # Exit with error code

    # Output JSON to stdout
    json.dump(data, sys.stdout, indent=2)

    if write_to_file:
        outpath = os.path.splitext(infile)[0] + ".json"
        try:
            with open(outpath, "w") as f:
                json.dump(data, f, indent=2)
            print(f"Successfully wrote JSON to {outpath}", file=sys.stderr)
        except Exception as e:
            print(f"Error writing to file {outpath}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
