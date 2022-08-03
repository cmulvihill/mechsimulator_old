import sys
from mechsimulator import runner

JOB_FILES = [
    # 'dme_sandia.xlsx',
    # 'dme_sandia_sens.xlsx',
    # 'dme_sandia_jsr_only.xlsx',
    # 'butane_jsr.xlsx'
    # 'madeup_h2o2_jsr.xlsx',
    # 'madeup_h2o2_st_idt.xlsx',
    # 'pham_n2o+o_check.xlsx'
    # 'n2o+o(1d)_tests.xlsx'
    # 'davidson_1990_sens.xlsx'
    # 'ro2_oh_idt_sens_dme.xlsx'
    # 'dme_sandia_515.xlsx',
    # 'c3h8_curran_tests.xlsx'
    # 'dme_sandia_513.xlsx',
    'dme_sandia_qooh_and_roo.xlsx',
    # 'c3h8_rcm_dames_2016.xlsx'
    # 'h2o2_test_rcm.xlsx'
]

runner.main.run_jobs(JOB_FILES)
