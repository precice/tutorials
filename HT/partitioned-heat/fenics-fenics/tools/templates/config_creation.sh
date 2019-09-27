{% for coupling_scheme in coupling_schemes -%}
{% for window_size in window_sizes -%}
{% for wr_left in wr_lefts -%}
{% for wr_right in wr_rights -%}
python3 create_waveform_config.py -exec {{executable}} -wr {{ wr_left }} {{ wr_right }} -dT {{ window_size }} -cpl {{ coupling_scheme }} --solver-tolerance {{solver_tolerance}} -qntol {{quasi_newton_tolerance}} -dd {{domain_decomposition}} -t {{time_dependence}} -mth {{method}} -wri {{waveform_interpolation_strategy}} -pp {{post_processing}} {{case_flag}} -T {{simulation_time}}
{% endfor -%}
{% endfor -%}
{% endfor -%}
{% endfor -%}
mv configuration.json experiments
