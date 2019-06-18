{% for window_size in window_sizes -%}
{% for wr_left in wr_lefts -%}
{% for wr_right in wr_rights -%}
cd experiments/WR{{wr_left}}{{wr_right}}/dT{{window_size}}/
./runall.sh
cd ../../..
{% endfor -%}
{% endfor -%}
{% endfor -%}
