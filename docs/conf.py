project = "Mass spectrometry-based Proteomics introduction"
copyright = "2025, DTU Biosustain, Informatics Platform, DSP"
author = "Henry Webel, Alberto Santos, Marco Reverenna"

extensions = [
    "myst_nb",
    "sphinx_new_tab_link",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/pandoc_ipynb/inputs/*",
    ".nox/*",
    "README.md",
    "**/.ipynb_checkpoints/*",
    "jupyter_execute",
    "conf.py",
]

html_theme = "pydata_sphinx_theme"

#  https://myst-nb.readthedocs.io/en/latest/computation/execute.html
nb_execution_mode = "off"

myst_enable_extensions = ["dollarmath", "amsmath"]

# Plotly support through require javascript library
# https://myst-nb.readthedocs.io/en/latest/render/interactive.html#plotly
html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js"
]

# https://myst-nb.readthedocs.io/en/latest/configuration.html
# Execution
nb_execution_raise_on_error = True
# Rendering
nb_merge_streams = True

# https://myst-nb.readthedocs.io/en/latest/authoring/custom-formats.html#write-custom-formats
nb_custom_formats = {".py": ["jupytext.reads", {"fmt": "py:percent"}]}

# https://myst-nb.readthedocs.io/en/latest/configuration.html#warning-suppression
suppress_warnings = ["mystnb.unknown_mime_type"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# See:
# https://github.com/executablebooks/MyST-NB/blob/master/docs/conf.py
# html_title = ""
html_theme = "sphinx_book_theme"
# html_logo = "_static/logo-wide.svg"
# html_favicon = "_static/logo-square.svg"
html_theme_options = {
    "github_url": "https://github.com/biosustain/dsp_course_proteomics_intro",
    "repository_url": "https://github.com/biosustain/dsp_course_proteomics_intro",
    # "repository_branch": "main",
    # "home_page_in_toc": True,
    # "path_to_docs": "docs",
    "show_navbar_depth": 1,
    # "use_edit_page_button": True,
    "use_repository_button": True,
    "use_download_button": True,
    "launch_buttons": {
        "colab_url": "https://colab.research.google.com"
        #     "binderhub_url": "https://mybinder.org",
        #     "notebook_interface": "jupyterlab",
    },
    "navigation_with_keys": False,
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
