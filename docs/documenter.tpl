{# This file based on https://github.com/marius311/CMBLensing.jl/blob/v0.10.1/docs/documenter.tpl
Under the MIT Expat license, © 2019–2023 Marius Millea #}

{% extends 'markdown/index.md.j2' %}


{% block stream %}
```output
{{ output.text | trim }}
```
{% endblock stream %}


{% block data_html scoped %}
```@raw html
{{ output.data['text/html'] }}
```
{% endblock data_html %}


{% block data_text scoped %}
```output
{{ output.data['text/plain'] | trim }}
```
{% endblock data_text %}