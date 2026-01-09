import glob
import re
import json


def main():
    generate_interfacing_snippets("python/*.py", "Python", "python/error_handling.py")
    generate_interfacing_snippets("rcpp/*.R", "R", "rcpp/error_handling.R")


def generate_interfacing_snippets(file_pattern, lang, *ignored):
    snippets = {}

    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored)):
            with open(file_name) as file:
                body = []
                description = None
                wrapped_func_name = None
                for line in file:
                    if line.startswith("#>"):
                        update_interfacing_snippets(snippets, body, wrapped_func_name, description, lang)
                        body = []
                        wrapped_func_name, description = get_interfacing_metadata(line)
                    else:
                        body.append(line)
                update_interfacing_snippets(snippets, body, wrapped_func_name, description, lang)

    print(json.dumps(snippets, indent=2))


def get_interfacing_metadata(initial_line):
    return re.search(r"^#> *(?P<wrapped_func_name>\w+):(?P<description>.*)", initial_line).groups()


def update_interfacing_snippets(snippets, body, wrapped_func_name, description, lang):
    if description is not None:
        func_name, args = get_interfacing_function_definition(body, lang)
        snippets.update(generate_definition_snippet(body, description, wrapped_func_name, lang))
        snippets.update(generate_call_snippet(func_name, wrapped_func_name, args, lang))


def generate_definition_snippet(body, description, wrapped_func_name, lang):
    return {
        f"{lang} setup for {wrapped_func_name}": {
            "prefix": f"tox:{lang[:2].lower()}_setup_{wrapped_func_name}",
            "body": "body",
            "description": description
        }
    }


def generate_call_snippet(func_name, wrapped_func_name, args, lang):
    match lang:
        case "Python":
            assignment_operator = "="
        case "R":
            assignment_operator = "<-"
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")

    arg_list_string = ", ".join(f"${{{i}:{arg}}}" for i, arg in enumerate(args, 1))
    return {
        f"{lang} call of {func_name}": {
            "prefix": f"tox:{lang[:2].lower()}_{func_name}",
            "body": f"result {assignment_operator} {func_name}({arg_list_string})",
            "description": f"Call of {lang} wrapper for {wrapped_func_name}"
        }
    }


def get_interfacing_function_definition(body, lang):
    """
    Extract (function_name, [args]) from a snippet body.
    `body` is a list of lines (strings).
    """

    match lang:
        case "Python":
            return _get_python_definition(body)
        case "R":
            return _get_r_definition(body)
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")


def _get_python_definition(body):
    # Remove inline comments and join lines
    cleaned = []
    for line in body:
        # strip inline comments
        line = re.sub(r"#.*", "", line)
        cleaned.append(line)

    text = "".join(cleaned)

    # def foo(a,b,c):
    m = re.search(r"def\s+([A-Za-z_]\w*)\s*\((.*?)\)\s*:", text, re.S)
    if not m:
        raise ValueError("No Python function definition found")

    func_name = m.group(1)
    arglist = m.group(2).strip()

    args = []
    if arglist:
        for a in arglist.split(","):
            a = a.strip()
            a = a.split("=")[0].strip()  # remove defaults
            a = a.lstrip("*")            # remove *args / **kwargs
            if a:
                args.append(a)

    return func_name, args


def _get_r_definition(body):
    cleaned = []
    for line in body:
        line = re.sub(r"#.*", "", line)
        cleaned.append(line)

    text = "".join(cleaned)

    # foo <- function(a,b,c)
    m = re.search(r"([A-Za-z_]\w*)\s*<-\s*function\s*\((.*?)\)", text, re.S)
    if not m:
        raise ValueError("No R function definition found")

    func_name = m.group(1)
    arglist = m.group(2).strip()

    args = []
    if arglist:
        for a in arglist.split(","):
            a = a.strip()
            a = a.split("=")[0].strip()
            if a:
                args.append(a)

    return func_name, args


def _get_fortran_definition(body):
    cleaned = []
    for line in body:
        # remove comments
        line = re.sub(r"!.*", "", line)
        # remove continuation markers
        line = line.replace("&", "")
        cleaned.append(line)

    text = "".join(cleaned).lower()

    # function foo(a,b,c)
    m = re.search(r"\b(function|subroutine)\s+([a-z_]\w*)\s*\((.*?)\)", text, re.S)
    if not m:
        raise ValueError("No Fortran function definition found")

    kind = m.group(1)
    func_name = m.group(2)
    arglist = m.group(3).replace("&", "").strip()

    args = []
    if arglist:
        for a in arglist.split(","):
            a = a.strip()
            if a:
                args.append(a)

    return func_name, args, kind


if __name__ == '__main__':
    main()
