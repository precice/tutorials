{% for window_size in window_sizes -%}
{% for wr_left in wr_lefts -%}
{% for wr_right in wr_rights -%}
{% for first_participant in first_participants -%}
cd experiments/WR{{wr_left}}{{wr_right}}/dT{{window_size}}/first_{{first_participant}}/
./runall.sh
cd ../../../..
{% endfor -%}
{% endfor -%}
{% endfor -%}
{% endfor -%}
