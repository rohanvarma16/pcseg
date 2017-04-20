FILE(REMOVE_RECURSE
  "CMakeFiles/hello.dir/./hello_generated_hello.cu.o"
  "hello.pdb"
  "hello"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/hello.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
