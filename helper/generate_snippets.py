import glob
import re
import json
from sys import stderr


def main():
    generate_interfacing_snippets("python/*.py", "Python", ["python/error_handling.py"])
    generate_interfacing_snippets("rcpp/*.R", "R", ["rcpp/error_handling.R"])
    generate_fortran_snippets("src/**/*.[fF]90", "Fortran_tox_snippets", ["src/f42_utils.F90", "src/config.F90", "src/safeguard.F90"])
    generate_fortran_snippets("src/f42_utils.F90", "Fortran_f42_snippets")


def generate_fortran_snippets(file_pattern, outfile_basename, ignored_files=[]):
    snippets = {}

    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                multiline = None
                module = None
                func_name = None
                args = None
                def_kind = None
                snippet_count = 0
                for line in file:
                    line = line.lower()
                    if re.search(r"end module", line.lower()) is not None:
                        break
                    if module is None:
                        module_def = re.search(r"^module\s+([a-z_1-9]+)", line)
                        if module_def is not None:
                            module, = module_def.groups()
                    else:
                        if not multiline:
                            func_def = re.search(r"\b(?P<def_kind>function|subroutine)\s+(?P<func_name>[a-z_0-9]\w*)\s*\((?P<args>[a-z_0-9,\s]*)(?P<continuation_char>&|\))?", line)
                            if func_def is not None:
                                def_kind, func_name, args, continuation_char = func_def.groups()
                                multiline = continuation_char == "&"
                        else:
                            more_args, continuation_char = re.search(r"\s*(?P<more_args>[a-z_0-9, ]*)(?P<continuation_char>&|\))?", line).groups()
                            args += more_args
                            multiline = continuation_char == "&"

                        if not multiline and multiline is not None:
                            call_str = "call " if def_kind == "subroutine" else "result = "
                            arg_list_string = arg_list_to_string(map(str.strip, args.split(",")))

                            snippets[f"Fortran call of {def_kind} {func_name}"] = {
                                "prefix": f"tox:{func_name}",
                                "body": f"{call_str} {func_name}({arg_list_string})",
                                "description": f"Call of Fortran {def_kind} {func_name} from module {module}"
                            }

                            func_name = None
                            args = None
                            def_kind = None
                            multiline = None
                            snippet_count += 1

                if module is None:
                    raise RuntimeError(f"No module found in file: {file_name}")
                if snippet_count == 0:
                    raise RuntimeError(f"No snippets generated for file: {file_name}")

    with open(f"snippets/{outfile_basename}.json", "w") as f:
        json.dump(snippets, f, indent=2)


def generate_interfacing_snippets(file_pattern, lang, ignored_files=[]):
    snippets = {}

    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                body = []
                description = None
                wrapped_func_name = None
                module = None
                indentation_state = "unset"
                for line_number, line in enumerate(file, 1):
                    if line.startswith("#>"):
                        update_interfacing_snippets(snippets, body, module, wrapped_func_name, description, lang)
                        body = []
                        module, wrapped_func_name, description = get_interfacing_metadata(line, file_name, line_number)
                        indentation = "unset"
                    elif module is not None:
                        # treat line as body (without line break)
                        body.append(line[:-1])

                        if len(line) > 1:
                            current_indentation = len(line) - len(line.lstrip(' '))

                            if current_indentation > 0:
                                if indentation_state == "unset":
                                    indentation_state = "indented"
                                elif indentation_state == "unindented":
                                    error_in_line("Declaring a new function without snippet header", file_name, line_number - 1)
                            elif indentation_state == "indented":
                                indentation_state = "unindented"

                update_interfacing_snippets(snippets, body, module, wrapped_func_name, description, lang)

    with open(f"snippets/{lang}_snippets.json", "w") as f:
        json.dump(snippets, f, indent=2)


def error_in_line(msg, file_name, line_number, post_msg=None):
    err_msg = f"{msg} in file '{file_name}', line {line_number}"
    if post_msg is not None:
        err_msg = f"{err_msg}: {post_msg}"
    raise RuntimeError(err_msg)


def get_interfacing_metadata(initial_line, file_name, line_number):
    metadata = re.search(r"^#> *(?P<module>[a-z_1-9A-Z]+):(?P<wrapped_func_name>[a-z_1-9A-Z]+):(?P<description>.*)", initial_line)
    if metadata is None:
        error_in_line("Invalid snippet header", file_name, line_number, "Should match pattern: '#> <fortran_module_name>:<fortran_subroutine_name>:<description>'")
    return metadata.groups()


def update_interfacing_snippets(snippets, body, module, wrapped_func_name, description, lang):
    if module is not None:
        func_name, args = get_interfacing_function_definition(body, lang)
        snippets.update(generate_definition_snippet(body, description, module, wrapped_func_name, lang))
        snippets.update(generate_interfacing_call_snippet(func_name, module, wrapped_func_name, args, lang))


def generate_definition_snippet(body, description, module, wrapped_func_name, lang):
    return {
        f"{lang} setup for {wrapped_func_name} from module '{module}'": {
            "prefix": f"tox:{lang[:2].lower()}_setup_{wrapped_func_name}",
            "body": body,
            "description": description
        }
    }


def arg_list_to_string(args):
    return ", ".join(f"${{{i}:{arg}}}" for i, arg in enumerate(args, 1))


def generate_interfacing_call_snippet(func_name, module, wrapped_func_name, args, lang):
    match lang:
        case "Python":
            assignment_operator = "="
        case "R":
            assignment_operator = "<-"
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")

    arg_list_string = arg_list_to_string(args)
    return {
        f"{lang} call of {func_name}": {
            "prefix": f"tox:{lang[:2].lower()}_{func_name}",
            "body": f"result {assignment_operator} {func_name}({arg_list_string})",
            "description": f"Call of {lang} wrapper for {wrapped_func_name} from module {module}"
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


def get_fortran_definition(body):
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
