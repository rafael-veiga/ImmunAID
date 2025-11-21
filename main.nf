#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CONSTRUCT_DATA {
  publishDir 'data_aux', mode: 'copy'

  input:
  path data_dir        
  path rscript         

  output:
  path 'data_raw.rds'     
  script:
  """
  mkdir -p data_aux
  Rscript ${rscript}
  """
}

process A1_PROCESS {
  publishDir 'data_aux', mode: 'copy'

  input:
  path data       
  path rscript        
  
  output:
  path 'data.rds'
  path 'data_nor.rds'
  path 'data_imp.rds'

  script:
  """
  Rscript ${rscript}
  """
}

process A2_PROCESS {

  input:
  path data_nor       
  path rscript        
  
  output:
  path 'mark_f'
  path 'mark_p'
  path 'data_nor_f.rds'
  path 'data_nor_p.rds'

  script:
  """
  mkdir mark_f
  mkdir mark_p
  Rscript ${rscript}
  """
}

process A3_PROCESS {

  input:
  path data_nor_f
  path data_nor_p      
  path script        
  
  output:
  path 'RF_auc1_f.csv'
  path 'RF_auc1_p.csv'
  path 'RF_auc2_f.csv'
  path 'RF_auc2_p.csv'
  path 'RF_imp_f.csv'
  path 'RF_imp_p.csv'
  path 'RF_imp2_f.csv'
  path 'RF_imp2_p.csv'

  script:
  """
  python ${script}
  """
}
process A4_PROCESS {

  input:
  path data_nor_f
  path data_nor_p
  path mark_f
  path mark_p   
  path script        
  
  output:
  path 'auc_n_f.csv'
  path 'auc_n_p.csv'
  path 'auc_curv_f.csv'
  path 'auc_curv_p.csv'
  
  script:
  """
  python ${script}
  """
}

process FIG1 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path auc_n_f
  path data
  path auc_curv_f
  path data_nor
  path rscript        
  
  output:
  path 'fig1.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG2 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path auc_n_f
  path data
  path auc_curv_f
  path data_nor
  path rscript        
  
  output:
  path 'fig2.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG3 {
  publishDir 'result', mode: 'copy'

  input:
  path RF_auc2_f
  path data_nor
  path rscript        
  
  output:
  path 'fig3.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG4 {
  publishDir 'result', mode: 'copy'

  input:
  path RF_imp2_f
  path data
  path rscript        
  
  output:
  path 'fig4.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG5 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_p
  path data
  path auc_n_p
  path auc_curv_p
  path data_nor
  path rscript        
  
  output:
  path 'fig5.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG6 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path mark_p
  path data_nor
  path rscript        
  
  output:
  path 'fig6.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG7 {
  publishDir 'result', mode: 'copy'

  input:
  path RF_auc2_p
  path data_nor
  path rscript        
  
  output:
  path 'fig7.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process FIG8 {
  publishDir 'result', mode: 'copy'

  input:
  path RF_imp2_p
  path data
  path rscript        
  
  output:
  path 'fig8.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S1 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path rscript        
  
  output:
  path 'figS1.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S2 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path auc_n_f
  path data
  path auc_curv_f
  path data_nor
  path rscript        
  
  output:
  path 'figS2.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S3 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path rscript        
  
  output:
  path 'figS3.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S4 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path auc_n_f
  path data
  path auc_curv_f
  path data_nor
  path rscript        
  
  output:
  path 'figS4.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S5 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path rscript        
  
  output:
  path 'figS5.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S6 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path auc_n_f
  path data
  path auc_curv_f
  path data_nor
  path rscript        
  
  output:
  path 'figS6.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S7 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_f
  path rscript        
  
  output:
  path 'figS7.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S8 {
  publishDir 'result', mode: 'copy'

  input:
  path RF_imp2_f
  path rscript        
  
  output:
  path 'figS8.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S9 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path mark_p
  path rscript        
  
  output:
  path 'figS9.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S10 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_p
  path data
  path auc_n_p
  path auc_curv_p
  path data_nor
  path rscript        
  
  output:
  path 'figS10.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S11 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path mark_p
  path rscript        
  
  output:
  path 'figS11.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S12 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_p
  path data
  path auc_n_p
  path auc_curv_p
  path data_nor
  path rscript        
  
  output:
  path 'figS12.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S13 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path mark_p
  path rscript        
  
  output:
  path 'figS13.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S14 {
  publishDir 'result', mode: 'copy'

  input:
  path mark_p
  path data
  path auc_n_p
  path auc_curv_p
  path data_nor
  path rscript        
  
  output:
  path 'figS14.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S15 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path mark_p
  path rscript        
  
  output:
  path 'figS15.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

process S16 {
  publishDir 'result', mode: 'copy'

  input:
  path data
  path RF_imp2_f
  path rscript        
  
  output:
  path 'figS16.pdf'
    
  script:
  """
  Rscript ${rscript}
  """
}

workflow {
  data_base_ch = channel.fromPath('data', type: 'dir')
  dara_raw_ch = CONSTRUCT_DATA( data_base_ch, channel.fromPath('scripts/Construct_datase_prot.R'))
  (data_ch, data_nor_ch, data_imp_ch) = A1_PROCESS(dara_raw_ch, channel.fromPath('scripts/A1_pre_process_data_analise.R'))
(mark_f_ch, mark_p_ch, data_nor_f_ch, data_nor_p_ch) = A2_PROCESS(data_nor_ch, channel.fromPath('scripts/A2_pre_process_data_analise.R'))
(RF_auc1_f_ch, RF_auc1_p_ch, RF_auc2_f_ch, RF_auc2_p_ch, RF_imp_f_ch, RF_imp_p_ch, RF_imp2_f_ch, RF_imp2_p_ch) = A3_PROCESS(data_nor_f_ch,data_nor_p_ch, channel.fromPath('scripts/analise4.py'))
(auc_n_f_ch, auc_n_p_ch, auc_curv_f_ch, auc_curv_p_ch) = A4_PROCESS(data_nor_f_ch, data_nor_p_ch, mark_f_ch, mark_p_ch, channel.fromPath('scripts/auc_curv.py'))
 FIG1(mark_f_ch, auc_n_f_ch, data_ch, auc_curv_f_ch, data_nor_ch, channel.fromPath('scripts/fig1.R'))
 FIG2(mark_f_ch, auc_n_f_ch, data_ch, auc_curv_f_ch, data_nor_ch, channel.fromPath('scripts/fig2.R'))
 FIG3(RF_auc2_f_ch, data_nor_ch, channel.fromPath('scripts/fig3.R'))
 FIG4(RF_imp2_f_ch, data_ch, channel.fromPath('scripts/fig4.R'))
 FIG5(mark_p_ch, data_ch, auc_n_p_ch, auc_curv_p_ch, data_nor_ch, channel.fromPath('scripts/fig5.R'))
 FIG6(data_ch, mark_p_ch, data_nor_ch, channel.fromPath('scripts/fig6.R'))
 FIG7(RF_auc2_p_ch, data_nor_ch, channel.fromPath('scripts/fig7.R'))
 FIG8(RF_imp2_p_ch, data_ch, channel.fromPath('scripts/fig8.R'))
 S1(mark_f_ch ,channel.fromPath('scripts/figS1.R'))
 S2(mark_f_ch, auc_n_f_ch, data_ch, auc_curv_f_ch, data_nor_ch, channel.fromPath('scripts/figS2.R'))
 S3(mark_f_ch ,channel.fromPath('scripts/figS3.R'))
 S4(mark_f_ch, auc_n_f_ch, data_ch, auc_curv_f_ch, data_nor_ch, channel.fromPath('scripts/figS4.R'))
 S5(mark_f_ch ,channel.fromPath('scripts/figS5.R'))
 S6(mark_f_ch, auc_n_f_ch, data_ch, auc_curv_f_ch, data_nor_ch, channel.fromPath('scripts/figS6.R'))
 S7(mark_f_ch ,channel.fromPath('scripts/figS7.R'))
 S8(RF_imp2_f_ch ,channel.fromPath('scripts/figS8.R'))

 S9(data_ch,mark_p_ch ,channel.fromPath('scripts/figS9.R'))
 S10(mark_p_ch, data_ch, auc_n_p_ch, auc_curv_p_ch, data_nor_ch, channel.fromPath('scripts/figS10.R'))
 S11(data_ch,mark_p_ch ,channel.fromPath('scripts/figS11.R'))
 S12(mark_p_ch, data_ch, auc_n_p_ch, auc_curv_p_ch, data_nor_ch, channel.fromPath('scripts/figS12.R'))
 S13(data_ch,mark_p_ch ,channel.fromPath('scripts/figS13.R'))
 S14(mark_p_ch, data_ch, auc_n_p_ch, auc_curv_p_ch, data_nor_ch, channel.fromPath('scripts/figS14.R'))
 S15(data_ch,mark_p_ch ,channel.fromPath('scripts/figS15.R'))
 S16(data_ch,RF_imp2_p_ch ,channel.fromPath('scripts/figS16.R'))

}
