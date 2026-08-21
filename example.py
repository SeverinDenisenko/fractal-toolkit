import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent / "python"))
    import fractaltoolkit as fk

    return fk, np


@app.cell
def _(fk, np):
    X = np.zeros(100)

    fk.generators.generate_white(X)
    fk.hurst.estimate_hurst_yw(X, 10)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
