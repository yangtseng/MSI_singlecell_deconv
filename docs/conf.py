# -- Sphinx configuration -------------------------------------------------
project   = "Single-cell MSI — PCMT cell typing"
author    = "Yang Tseng"
extensions = ["sphinx.ext.mathjax"]      # math directives; no third-party extension needed
master_doc = "index"
exclude_patterns = ["_build"]

# Use the ReadTheDocs theme if installed, else fall back to a built-in one.
try:
    import sphinx_rtd_theme            # noqa: F401
    html_theme = "sphinx_rtd_theme"
except ImportError:
    html_theme = "alabaster"           # ships with Sphinx
