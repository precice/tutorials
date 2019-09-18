{% for coupling_scheme in coupling_schemes -%}
{% for window_size in window_sizes -%}
{% for wr_left in wr_lefts -%}
{% for wr_right in wr_rights -%}
python3 create_waveform_config.py -exec {{executable}} -wr {{ wr_left }} {{ wr_right }} -dT {{ window_size }} -cpl {{ coupling_scheme }} --gamma {{gamma}} --solver-tolerance {{solver_tolerance}} -qntol {{quasi_newton_tolerance}} -dd {{domain_decomposition}} {{plain_subcycling}} -t {{time_dependence}} -mth {{method}} -wri {{waveform_interpolation_strategy}} -pp {{post_processing}}
{% endfor -%}
{% endfor -%}
{% endfor -%}
{% endfor -%}
