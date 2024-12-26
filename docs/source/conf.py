# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Cradar'
copyright = '2022, sfranke'
author = 'sfranke'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    "nbsphinx",
]



intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

extensions.append("sphinx_wagtail_theme")
html_theme = 'sphinx_wagtail_theme'

pygments_style = 'sphinx'

# -- Options for EPUB output
epub_show_urls = 'footnote'

html_sidebars = {
    "**": ["sidebar-nav-bs", "sidebar-ethical-ads"]
}
