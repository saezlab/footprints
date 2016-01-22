.io = import('io')

files = list.files(module_file(), "\\.yaml$")

for (.file in files)
    assign(sub("\\.yaml", "", .file), .io$read_yaml(.file))
