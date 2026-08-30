genrule(
    name = "main",
    srcs = [
        "//src/RFL:rfl_wheel",
        "main.py",
    ],
    outs = ["run_main.sh"],
    cmd = [
        "echo '#!/bin/bash' > $OUTS",
        "echo 'python3 -m venv .main_venv' >> $OUTS",
        "echo 'source .main_venv/bin/activate' >> $OUTS",
        "echo 'pip install --quiet numpy' >> $OUTS",
        "echo 'pip install --quiet --no-index --find-links $(out_location //src/RFL:rfl_wheel) rfl' >> $OUTS",
        "echo 'python src/RFL/examples/python/main.py' >> $OUTS",
        "chmod +x $OUTS",
    ],
    binary = True,
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
        "echo 'pip install --quiet --no-index --find-links $(out_location //src/RFL:rfl_wheel) rfl' >> $OUTS",
        "echo 'echo \"Launching Jupyter Notebook...\"' >> $OUTS",
        "echo 'jupyter notebook src/RFL/examples/python/rfl_playground.ipynb' >> $OUTS",
        "chmod +x $OUTS",
    ],
    binary = True,
    visibility = ["PUBLIC"],
)
