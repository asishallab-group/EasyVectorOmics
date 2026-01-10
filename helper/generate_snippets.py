import glob
import re
import json
from sys import stderr

def main():
    fortran_wrappers, fortran_modules = generate_fortran_snippets("src/**/*.[fF]90", ["src/config.F90", "src/safeguard.F90"])
    generate_interfacing_snippets("python/*.py", "Python", fortran_modules, fortran_wrappers={*fortran_wrappers})
    generate_interfacing_snippets("rcpp/*.R", "R", fortran_modules, fortran_wrappers={*fortran_wrappers})


def get_initial_snippet_dicts():
    return {"_kind": "tox"}, {"_kind": "f42"}


def generate_fortran_snippets(file_pattern, ignored_files=[]):
    tox_snippets, f42_snippets = get_initial_snippet_dicts()
    fortran_wrappers = set()
    fortran_modules = set()

    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                multiline = None
                module = None
                module_snippets = None
                func_name = None
                args = None
                def_kind = None
                snippet_count = 0
                is_wrapper_section = False
                for line in file:
                    line = line.lower()
                    if re.search(r"end module", line.lower()) is not None:
                        is_wrapper_section = True
                    if module is None:
                        module_def = re.search(r"^module\s+([a-z_0-9]+)", line)
                        if module_def is not None:
                            module, = module_def.groups()
                            fortran_modules.add(module.lower())
                            module_snippets = get_module_snippets(module, f42_snippets, tox_snippets)
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
                            if is_wrapper_section and func_name[-1].lower() == "c":
                                fortran_wrappers.add(func_name.lower())
                            else:
                                call_str = "call" if def_kind == "subroutine" else "${{0:result}} ="
                                arg_list_string = arg_list_to_string(map(str.strip, args.split(",")))

                                module_snippets[f"Fortran Call of {def_kind} '{func_name}'"] = {
                                    "prefix": f"{module_snippets["_kind"]}:{func_name}",
                                    "body": f"{call_str} {func_name}({arg_list_string})",
                                    "description": f"Call of Fortran {def_kind} '{func_name}' from module '{module}'"
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

    write_snippets(f42_snippets, tox_snippets, "Fortran")
    return fortran_wrappers, fortran_modules


def get_module_snippets(module_name, f42_snippets, tox_snippets):
    if module_name.startswith("f42"):
        return f42_snippets
    if module_name.startswith("tox"):
        return tox_snippets
    else:
        raise RuntimeError(f"module name '{module_name}' does not start with 'tox' or 'f42'")


def generate_interfacing_snippets(file_pattern, lang, fortran_modules, ignored_files=[], fortran_wrappers=set()):
    tox_snippets, f42_snippets = get_initial_snippet_dicts()
    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                body = []
                description = None
                wrapped_func_name = None
                module = None
                module_snippets = None
                defined_function = False
                is_void = True
                for line_number, line in enumerate(file, 1):
                    if line.startswith("#>"):
                        if line.startswith("#>skip snippets"):
                            break
                        update_interfacing_snippets(module_snippets, body, module, wrapped_func_name, description, is_void, lang, fortran_wrappers, file_name, line_number)
                        is_void = True
                        body = []
                        module, wrapped_func_name, description = get_interfacing_metadata(line, file_name, line_number)
                        defined_function = False
                        if not module.endswith(":helper"):
                            if module.lower() not in fortran_modules:
                                raise error_in_line(f"Unknown module '{module}'", file_name, line_number)
                            if wrapped_func_name.lower() not in fortran_wrappers:
                                raise error_in_line(f"Unknown interfacing subroutine '{wrapped_func_name}' for module '{module}'", file_name, line_number, "(Could also be duplicate)")
                        elif wrapped_func_name is not None:
                            defined_function = True
                        module_snippets = get_module_snippets(module, f42_snippets, tox_snippets)
                    elif module is not None:
                        if re.search(r'(<-\s*function|^def)', line) is not None:
                            if defined_function:
                                error_in_line("Declaring a new function without snippet header (or just missing indentation)", file_name, line_number)
                            defined_function = True

                        # treat line as body (without line break)
                        body.append(line[:-1])

                        if re.search(r"^\s*\breturn\b", line) is not None:
                            is_void = False

                update_interfacing_snippets(module_snippets, body, module, wrapped_func_name, description, is_void, lang, fortran_wrappers, file_name, line_number)

    if lang == "Python" and len(fortran_wrappers) > 0:
        raise RuntimeError(f"Missing Python wrappers for: {", ".join(fortran_wrappers)}")

    write_snippets(f42_snippets, tox_snippets, lang)


def write_snippets(f42_snippets, tox_snippets, lang):
    for snippets in ["tox_snippets", "f42_snippets"]:
        with open(f"snippets/{lang}_{snippets}.json", "w") as f:
            snippets_dict = locals()[snippets]
            kind = snippets_dict["_kind"]
            del snippets_dict["_kind"]
            json.dump(snippets_dict, f, indent=2)
            snippets_dict["_kind"] = kind


def error_in_line(msg, file_name, line_number, post_msg=None):
    err_msg = f"{msg} in file '{file_name}', line {line_number}"
    if post_msg is not None:
        err_msg = f"{err_msg}: {post_msg}"
    raise RuntimeError(err_msg)


def get_interfacing_metadata(initial_line, file_name, line_number):
    metadata = re.search(r"^#> *(?P<module>[a-z_0-9A-Z]+):(?P<wrapped_func_name>[a-z_0-9A-Z]+):\s*(?P<description>.*)\s*", initial_line)
    if metadata is not None:
        return metadata.group("module"), metadata.group("wrapped_func_name"), metadata.group("description")

    metadata = re.search(r"^#> (?P<kind>f42|tox)_helper(-(?P<name>[a-zA-Z0-9_]+))?:\s*(?P<description>.*)\s*", initial_line)
    if metadata is not None:
        return f"{metadata.group("kind")}:helper", metadata.group("name"), metadata.group("description")

    error_in_line("Invalid snippet header", file_name, line_number, "Should match pattern: '#> <fortran_module_name>:<fortran_subroutine_name>:<description>' OR '#> (tox|f42)_helper[-<name>]:<description>' ")


def update_interfacing_snippets(snippets, body, module, wrapped_func_name, description, is_void, lang, fortran_wrappers, file_name, line_number):
    if module is not None:
        if module.endswith(":helper") and wrapped_func_name is not None:
            remove_unnecessary_lines(body)
            snippets[f"{lang} Helper to {description}"] = {
                "prefix": f"{snippets["_kind"]}:{lang[:2].lower()}_helper_{wrapped_func_name}",
                "body": body,
                "description": description
            }
        else:
            func_name, args = get_interfacing_function_definition(body, lang)
            if module.endswith(":helper"):
                wrapped_func_name = func_name
            elif lang == "Python":
                fortran_wrappers.discard(wrapped_func_name.lower())

            snippets.update(generate_definition_snippet(body, description, module, wrapped_func_name, lang, snippets["_kind"]))
            snippets.update(generate_interfacing_call_snippet(func_name, module, wrapped_func_name, args, is_void, lang, snippets["_kind"]))


def remove_unnecessary_lines(body):
    # remove last unnecessary lines, either empty (or spaces) or comments
    while re.search(r'^\s*(#|"""|$)', (last_line := body.pop())) is not None:
        pass

    body.append(last_line)
    body.append("")


def generate_definition_snippet(body, description, module, wrapped_func_name, lang, snippet_kind):
    if module.endswith(":helper"):
        key = f"{lang} Helper function '{wrapped_func_name}'"
        prefix = f"{snippet_kind}:{lang[:2].lower()}_helper_{wrapped_func_name}"
    else:
        key = f"{lang} setup for '{wrapped_func_name}' from module '{module}'"
        prefix = f"{snippet_kind}:{lang[:2].lower()}_setup_{wrapped_func_name}"

    remove_unnecessary_lines(body)

    return {
        key: {
            "prefix": prefix,
            "body": body,
            "description": description
        }
    }


def arg_list_to_string(args, tab_index=1):
    return ", ".join(f"${{{i}:{arg}}}" for i, arg in enumerate(args, tab_index))


def generate_interfacing_call_snippet(func_name, module, wrapped_func_name, args, is_void, lang, snippet_kind):
    match lang:
        case "Python":
            assignment_operator = "="
        case "R":
            assignment_operator = "<-"
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")

    if module.endswith(":helper"):
        description = f"Call of {lang} helper '{wrapped_func_name}'"
    else:
        description = f"Call of {lang} wrapper for '{wrapped_func_name}' from module '{module}'"

    if not is_void:
        body = f"${{1:result}} {assignment_operator} "
        arg_list_string = arg_list_to_string(args, 2)
    else:
        body = ""
        arg_list_string = arg_list_to_string(args, 1)

    body += f"{func_name}({arg_list_string})"

    return {
        f"{lang} call of '{func_name}'": {
            "prefix": f"{snippet_kind}:{lang[:2].lower()}_{func_name}",
            "body": body,
            "description": description
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
    m = re.search(r"def\s+([A-Za-z_]\w*)\s*\((.*?)\)\s*(->\s*\w*\s*)?:", text, re.S)
    if not m:
        raise ValueError("No Python function definition found in:\n\n", text)

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


if __name__ == '__main__':
    main()
