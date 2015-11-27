

export GLIBCXX_FORCE_NEW

valgrind -v --num-callers=20 --leak-resolution=high --leak-check=yes --log-file="logfile.txt" --read-var-info=yes build/run_tests
# --leak-check=full --show-leak-kinds=all
