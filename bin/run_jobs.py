from mechsimulator import runner

JOB_FILES = [
    # 'dme_sandia.xlsx',
    'dme_sandia_sens.xlsx',
    # 'dme_sandia_jsr_only.xlsx',
    # 'butane_jsr.xlsx'
    # 'madeup_h2o2_jsr.xlsx',
    # 'pham_n2o+o_check.xlsx'
    # 'n2o+o(1d)_tests.xlsx'
    # 'davidson_1990_sens.xlsx'
    # 'ro2_oh_idt_sens_dme.xlsx'
    # 'dme_sandia_515.xlsx',
]

runner.main.run_jobs(JOB_FILES)
