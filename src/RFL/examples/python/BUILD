python_binary(
    name = "main",
    main = "main.py",
    deps = [
        "//src/RFL:rfl_wheel",
        "//third_party/python:numpy",
    ],
    visibility = ["PUBLIC"],
)

genrule(
    name = "jupyter",
    srcs = ["//src/RFL:rfl_wheel"],
    outs = ["run_jupyter.sh"],
    cmd = [
        "echo '#!/bin/bash' > $OUTS",
        "echo 'echo \"Setting up Jupyter virtual environment...\"' >> $OUTS",
        "echo 'python3 -m venv .jupyter_venv' >> $OUTS",
        "echo 'source .jupyter_venv/bin/activate' >> $OUTS",
        "echo 'pip install --quiet jupyter matplotlib numpy' >> $OUTS",
        "echo \"pip install --quiet $SRCS\" >> $OUTS",
        "echo 'echo \"Launching Jupyter Notebook...\"' >> $OUTS",
        "echo 'jupyter notebook src/RFL/examples/python/rfl_playground.ipynb' >> $OUTS",
        "chmod +x $OUTS"
    ],
    binary = True,
    visibility = ["PUBLIC"],
)
