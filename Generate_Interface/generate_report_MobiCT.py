#######################################################
#                  FUNCTIONS                          #
#######################################################
import re
import json
import pandas as pd
import numpy as np
import os
import codecs
import argparse
import subprocess

def whiteAnnot(p1l,gene,path):
    if p1l == "":
        return ""
    database = pd.read_csv(path,sep="\t",encoding='latin-1')
    database = database.loc[database["Gene"]==gene]
    if sum(p1l in i for i in database['Alterations'])>0:
        res = 'Yes'#np.array(database[[p1l in i for i in database['Alterations']]]['Level'])[0]
        return res
    else:
        return ""

def command(cmd):
    subprocess.check_call(cmd, shell=True)

def threeToOne(string):
    string = string.replace('Ala','A')
    string = string.replace('Asx','B')
    string = string.replace('Cys','C')
    string = string.replace('Asp','D')
    string = string.replace('Glu','E')
    string = string.replace('Phe','F')
    string = string.replace('Gly','G')
    string = string.replace('His','H')
    string = string.replace('Ile','I')
    string = string.replace('Lys','K')
    string = string.replace('Leu','L')
    string = string.replace('Met','M')
    string = string.replace('Asn','N')
    string = string.replace('Pro','P')
    string = string.replace('Gln','Q')
    string = string.replace('Arg','R')
    string = string.replace('Ser','S')
    string = string.replace('Thr','T')
    string = string.replace('Sec','U')
    string = string.replace('Val','V')
    string = string.replace('Trp','W')
    string = string.replace('Xaa','X')
    string = string.replace('Tyr','Y')
    string = string.replace('Glx','Z')
    string = string.replace('Ter','*')
    return(string)



def getTemplate():
    return """<!DOCTYPE html>
    <html>
    <head>
      <title>##sample_name##</title>
      <meta charset="UTF-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1" />
      <link rel="stylesheet" href="../ressources/lib///w3.css" />
      <link rel="stylesheet" href="../ressources/lib/fontawesome-free-5.3.1-web/css/all.css" />
      <link rel="stylesheet" href="../ressources/lib///datatables/datatables.min.css" />
      <link rel="stylesheet" href="../ressources/lib///interface.css" />
    </head>
    <body>
      <!-- Top container -->
      <div class="w3-bar w3-top w3-black w3-large" style="z-index:4">
        <button class="w3-bar-item w3-button w3-hide-large w3-hover-none w3-hover-text-light-grey" onclick="w3_open();"><i class="fa fa-bars"></i> Menu</button>
        <a href="run.html" class="w3-left w3-bar-item w3-button w3-teal">
          <i class="fa fa-home"></i>
        </a>
        <span class="w3-display-middle w3-bar-item w3-hide-small"><i class="fa fa-user"></i> ##sample_name##</span>
        <a href="run.html" class="w3-right w3-bar-item w3-button w3-red w3-hide-medium w3-hide-small"><i class="fa fa-digital-tachograph"></i>Back to Run</a>
        <a href="run.html" class="w3-right w3-bar-item w3-button w3-red w3-hide-large">
          <i class="fa fa-digital-tachograph"></i>
        </a>
      </div>
      <div id="body" class="w3-main" style="margin-left:30px;margin-top:43px;">
        <header id="sample" class="w3-container" style="padding-top:43px">
          <h1><i class="fa fa-user"></i> ##sample_name##</h1>
        </header>
        <body>
        <div id="Summary" class="w3-container" style="">
        <h5>
        <b><i class="fa fa-tachometer-alt"></i> Summary</b>
        </h5>
        <table id="Summary" class="w3-table w3-card" style="width:35%;margin-left:20px">
        <tr>
          <td>Limit of detection (p<0.01)</td>
          <td>##LoD##</td>
        </tr>
        <tr>
          <td>Total Deduplicated Reads </td>
          <td>##Total##</td>
        </tr>
        <tr>
          <td>Deduplication Rate </td>
          <td>##Dedup##</td>
        </tr>
        <tr>
          <td>Panel coverage </td>
          <td>##cov##</td>
        </tr>
        <tr>
          <td>Number of variants</td>
          <td>##number##</td>
        </tr>
      </table>
      </body>
      </div>
        <div id="Variants" class="w3-container" style="padding-top:43px;">
          <h5>
            <b><i class="fa fa-dna"></i> Variants</b>
          </h5>
          <div id="" class="w3-container" style="padding:0px;">
            <div id="sssamples" class="w3-padding" style="padding-top:430px;">
              <table id="T_all_variants" class="w3-table w3-hoverable" style="">
                <thead class="" style="">
                  <tr class="" style="">
                    <th class="" style="">Sample name</th>
                    <th class="" style="">Gene</th>
                    <th class="" style="">DNA position</th>
                    <th class="" style="">Protein position (1 letter)</th>
                    <th class="" style="">Protein position (3 letters)</th>
                    <th class="" style="">Localization</th>
                    <th class="" style="">VAF</th>
                    <th class="" style="">Consequences</th>
                    <th class="" style="">Filters</th>
                    <th class="" style="">Clinvar</th>
                    <th class="" style="">Whitelist</th>
                    <th class="" style="">VD</th>
                    <th class="" style="">DP</th>
                    <th class="" style="">chrom</th>
                    <th class="" style="">pos</th>
                    <th class="" style="">ref</th>
                    <th class="" style="">alt</th>
                    <th class="" style="">NM</th>
                  </tr>
                </thead>
                <tbody class="" style="">
                  <tr class="##color##" style="height:40px">
                    <td class="" style="">##sample_name##</td>
                    <td class="" style="">##gene##</td>
                    <td class="" style="">##dna_pos##</td>
                    <td class="" style="">##prot_pos1l##</td>
                    <td class="" style="">##prot_pos##</td>
                    <td class="" style="">##localization##</td>
                    <td class="" style="">##vaf##</td>
                    <td class="" style="">##csq##</td>
                    <td class="" style="">##filters##</td>
                    <td class="" style="">##clinvar##</td>
                    <td class="" style="">##white##</td>
                    <td class="" style="">##vd##</td>
                    <td class="" style="">##dp##</td>
                    <td class="" style="">##chrom##</td>
                    <td class="" style="">##pos##</td>
                    <td class="" style="">##ref##</td>
                    <td class="" style="">##alt##</td>
                    <td class="" style="">##nm##</td>
                  </tr>
                <tr class="w3-pale-orange w3-hover-orange na hidden" style="height:40px">
                  <td class="" style="">##sample_name##</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">0</td>
                  <td class="" style="">0</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                  <td class="" style="">NA</td>
                </tr>
                <tr class="w3-pale-orange w3-hover-orange wt hidden" style="height:40px">
                  <td class="" style="">##sample_name##</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">0</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                  <td class="" style="">WT</td>
                </tr>
                </tbody>
                <tfoot class="" style="">
                  <tr class="" style="">
                    <th class="" style="">Sample name</th>
                    <th class="" style="">Gene</th>
                    <th class="" style="">DNA position</th>
                    <th class="" style="">Protein position (1 letter)</th>
                    <th class="" style="">Protein position (3 letters)</th>
                    <th class="" style="">Localization</th>
                    <th class="" style="">VAF</th>
                    <th class="" style="">Consequences</th>
                    <th class="" style="">Filters</th>
                    <th class="" style="">Clinvar</th>
                    <th class="" style="">Whitelist</th>
                    <th class="" style="">VD</th>
                    <th class="" style="">DP</th>
                    <th class="" style="">chrom</th>
                    <th class="" style="">pos</th>
                    <th class="" style="">ref</th>
                    <th class="" style="">alt</th>
                    <th class="" style="">NM</th>
                  </tr>
                </tfoot>
              </table>
            </div>
          </div>
        </div>
      </div>
      <script type="text/javascript" src="../ressources/lib///interface.js"></script>
      <script type="text/javascript" src="../ressources/lib///datatables/datatables.min.js"></script>
      <script type="text/javascript" src="../ressources/lib///node_modules/highcharts/highstock.js"></script>
      <script type="text/javascript" src="../ressources/lib///d3js/d3.v3.min.js"></script>
      <script type="text/javascript" src="../ressources/lib///node_modules/highcharts/highcharts-more.js"></script>
      <script type="text/javascript" src="../ressources/lib///node_modules/highcharts/modules/exporting.js"></script>
      <script type="text/javascript" src="../ressources/lib///node_modules/highcharts/modules/export-data.js"></script>
      <!-- table filters -->
      <script>
              $(document).ready(function() {
                  $('#T_filters tfoot th').each( function () {
                      var title = $(this).text();
                      $(this).html( '<input type="text" placeholder="Search '+title+'" style="width:100%" />' );
                  } );

                  var table = $('#T_filters').DataTable({
                  scrollY: "400px",
                      scrollCollapse: false,
                      scrollX:false,
                      paging:false,
                      "dom": '<"top">t<"bottom"i><"clear">',
                      "order": [[ 0, "asc" ]]
                  });

                  table.columns().every( function () {
                      var that = this;

                      $( 'input', this.footer() ).on( 'keyup change', function () {
                          if ( that.search() !== this.value ) {
                              that
                                  .search( this.value )
                                  .draw();
                          }
                      });
                  });
              });
          </script>
      <!-- table amplicons -->
      <script>
              $(document).ready(function() {
                  $('#T_amplicons tfoot th').each( function () {
                      var title = $(this).text();
                      $(this).html( '<input type="text" placeholder="Search '+title+'" style="width:100%" />' );
                  } );

                  var table = $('#T_amplicons').DataTable({
                  scrollY: '400px',
                      scrollCollapse: false,
                      scrollX:false,
                      paging:false,
                      "dom": '<"top">t<"bottom"i><"clear">',
                  });

                  table.columns().every( function () {
                      var that = this;

                      $( 'input', this.footer() ).on( 'keyup change', function () {
                          if ( that.search() !== this.value ) {
                              that
                                  .search( this.value )
                                  .draw();
                          }
                      });
                  });
              });
          </script>
      <!-- table amplicons -->
      <script>
              $(document).ready(function() {
                  $('#T_all_variants tfoot th').each( function () {
                      var title = $(this).text();
                      $(this).html( '<input type="text" placeholder="Search '+title+'" style="width:100%" />' );
                  } );

                  var table = $('#T_all_variants').DataTable({
                      scrollY: "400px",
                      scrollX: true,
                      scrollCollapse: true,
                      paging: false,
                      "order": [[ 10, "desc" ],[ 6, "desc" ]],
                      select: {
                          style: 'multi',
                      },
                      "dom": '<"top"B>t<"bottom"><"clear">',
                      buttons: {
                          dom: {
                              button: {
                                  className: 'w3-button'
                              }
                          },
                          buttons: [
                              {
                                  extend: 'copy',
                                  title: '',
                                  text: 'Copy (NA)',
                                  exportOptions: {
                                      rows: '.na',
                                  },
                                  className: 'w3-orange',
                                  header: false
                              },
                              {
                                  extend: 'copy',
                                  title: '',
                                  text: 'Copy (WT)',
                                  exportOptions: {
                                      rows: '.wt',
                                  },
                                  className: 'w3-orange',
                                  header: false
                              },
                              {
                                  extend: 'copyHtml5',
                                  title: '',
                                  text: 'Copy (selected)',
                                  className: 'w3-blue',
                                  exportOptions: {
                                      modifier: {
                                          selected: true
                                      }
                                  },
                          header: false
                              },
                              {
                                  extend: 'excelHtml5',
                                  title: 'variants_selected_##sample_name##',
                                  text: 'Excel (selected)',
                                  className: 'w3-blue',
                                  exportOptions: {
                                      modifier: {
                                          selected: true
                                      }
                                  }
                              },
                              {
                                  extend: 'csvHtml5',
                                  title: 'variants_selected_##sample_name##',
                                  text: 'CSV (selected)',
                                  className: 'w3-blue',
                                  exportOptions: {
                                      modifier: {
                                          selected: true
                                      }
                                  }
                              },
                              {
                                  extend: 'copyHtml5',
                                  title: '',
                                  text: 'Copy (all)',
                                  className: 'w3-teal',
                          header: false,
                                  exportOptions: {
                                      rows: '.var',
                                      modifier: {
                                          selected: null
                                      }
                                  },
                              },
                              {
                                  extend: 'excelHtml5',
                                  title: 'variants_##sample_name##',
                                  text: 'Excel (all)',
                                  className: 'w3-teal',
                                  exportOptions: {
                                      rows: '.var',
                                      modifier: {
                                          selected: null
                                      }
                                  },
                              },
                              {
                                  extend: 'csvHtml5',
                                  title: 'variants_##sample_name##',
                                  text: 'CSV (all)',
                                  className: 'w3-teal',
                                  exportOptions: {
                                      rows: '.var',
                                      modifier: {
                                          selected: null
                                      }
                                  },
                              }
                          ]
                      }
                  });
                  table.columns().every( function () {
                      var that = this;

                      $( 'input', this.footer() ).on( 'keyup change', function () {
                          if ( that.search() !== this.value ) {
                              that
                                  .search( this.value )
                                  .draw();
                          }
                      });
                  });
              });

          </script>
    <script>
                   $('.hidden').hide()
                  </script>
    </body>
    </html>
"""


#######################################################
#                       MAIN                          #
#######################################################
parser = argparse.ArgumentParser(description='This function is used to create a Adivar-like report for Twist cfDNA runs')
parser.add_argument('-v', '--variants', type=str, help='path to the directory containing the vcf files')
parser.add_argument('-s', '--stats', type=str, help='path to to the directory containing the metrics')
parser.add_argument('-k', '--kpath', type=str, help='path to to the cov_vaf_probs csv file')
parser.add_argument('-w', '--whitelist', type=str, help='path to to the whitelist txt file')
parser.add_argument('-b', '--bed', type=str, help='path to to the bedfile')
parser.add_argument('-t', '--template', type=str, help='path to to the template.html file')
parser.add_argument('-o', '--output', type=str, help='output path')
args = parser.parse_args()

k_path = args.kpath
whitelist_file = args.whitelist

to_keep = pd.read_csv(args.bed,sep="\t",header=None)
to_keep = np.unique(to_keep[3]).tolist()
print(to_keep)
kept_csq = ['TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant', 'regulatory_region_ablation', 'regulatory_region_amplification', 'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'splice_region_variant']

k3 = pd.read_csv(k_path, sep="\t")
command = "mkdir -p "+args.output
subprocess.run(command, shell=True)


path = args.variants
samples=[]
for i in os.listdir(path):
    if i.endswith(".vcf"):
        samples.append(i.split("_vep.vcf")[0])

samples = np.unique(samples)

x="""                  <tr class="##color##" style="height:40px">
                    <td class="" style="">##sample_name##</td>
                    <td class="" style="">##gene##</td>
                    <td class="" style="">##dna_pos##</td>
                    <td class="" style="">##prot_pos1l##</td>
                    <td class="" style="">##prot_pos##</td>
                    <td class="" style="">##localization##</td>
                    <td class="" style="">##vaf##</td>
                    <td class="" style="">##csq##</td>
                    <td class="" style="">##filters##</td>
                    <td class="" style="">##clinvar##</td>
                    <td class="" style="">##white##</td>
                    <td class="" style="">##vd##</td>
                    <td class="" style="">##dp##</td>
                    <td class="" style="">##chrom##</td>
                    <td class="" style="">##pos##</td>
                    <td class="" style="">##ref##</td>
                    <td class="" style="">##alt##</td>
                    <td class="" style="">##nm##</td>
                  </tr>"""

for sample_name in samples:
    print("Doing " + sample_name + ":")

    # Set the output report
    output_report = args.output+"/"+sample_name+".html"

    # Gather statistics on the sample
    picard = pd.read_csv(args.stats+"/"+sample_name+"_output_hs_metrics1.txt",sep='\t',engine='python',skiprows=lambda x: x not in [6,7])
    metrics = pd.read_csv(args.stats+"/"+sample_name+"_reportAfter/multiqc_data/multiqc_picard_HsMetrics.txt",sep="\t")

    m = round(metrics['MIN_TARGET_COVERAGE'][0])
    M = round(metrics['MAX_TARGET_COVERAGE'][0])
    ave = round(metrics['MEAN_TARGET_COVERAGE'][0])
    dedup_rate = round(ave*100/picard['MEAN_TARGET_COVERAGE'][0],1)
    N = round(metrics['ON_TARGET_BASES'][0])
    if ave<200:
        LoD = 50
    else:
        LoD = min(k3.loc[[i<ave for i in k3['p99']]]['vaf'])

    # Get variants annotated with vep
    subprocess.run("grep -v \"^##\" "+args.variants+"/"+sample_name+"_vep.vcf > "+args.variants+"/"+sample_name+".tsv", shell=True, check=True)
    variant_data = pd.read_csv(args.variants+"/"+sample_name+".tsv",sep="\t")

    # Build output report
    html = getTemplate()

    ## header
    total = str(round(N/1000000,2))+"M reads"
    cov2 = "Min = "+str(m)+" | Max = "+str(M)+" | Average = "+str(round(ave))
    nvar = str(len(variant_data))

    html = html.replace("##sample_name##",sample_name,3)
    html = html.replace("##LoD##",str(LoD)+"%",1)
    html = html.replace("##Total##",total,1)
    html = html.replace("##Dedup##",str(dedup_rate)+'%',1)
    html = html.replace("##cov##",cov2,1)

    q=0
    for i in variant_data.index:
        tmp = variant_data.loc[i]

        whitelist = ""
        prot_pos = ""
        nm = ""
        p3l = ""
        p1l = ""

        filt = []
        filters = tmp[6]
        if filters == 'f0.02':
            filters='PASS'
        elif 'f0.02' in filters:
            filters = filters.split(";")[1]
        if filters == 'p8':
            filters = 'PASS'

        chrom = tmp[0]
        pos = tmp[1]
        ref = tmp[3]
        alt = tmp[4]
        filt = tmp[6]
        details = tmp[7].split(";")
        #print("------------------------------------")
        #print(i)
        #print(details)
        #if i==295:
        #    print(details)

        dp = np.array(details)[["DP=" in i for i in details]][0].split("DP=")[1]
        vd = np.array(details)[["VD=" in i for i in details]][0].split("VD=")[1]
        af = np.array(details)[["AF=" in i for i in details]][0].split("AF=")[1]

        if sum(["CSQ=" in i for i in details])>0:
            var = np.array(details)[["CSQ=" in i for i in details]][0].split(",")

            good = "x"
            for k in to_keep:
                for l in var:
                    if k in l:
                        good = k
            if good != "x":
                var = np.array(var)[[good in i for i in var]][0].split("|")
                gene = var[3]

                if ":" in var[10]:
                    c = var[10].split(':')[1]
                else:
                    c=''
                clinvar = var[32].replace('&', ', ')
                csq = var[1]
                keep_csq = 0
                for k in csq.split('&'):
                    if k in kept_csq:
                        keep_csq = keep_csq + 1
                csq = csq.replace('&',', ')

                if len(var[11])>0:
                    p3l = var[11].split(':')[1].split('.')[1]
                    p3l = p3l.replace("%3D","=")
                    p1l = threeToOne(p3l)
                if var[30]!='':
                    if float(var[30])>0.01:
                        filters = 'popAF'

                if var[8]=='':
                    loc = "Intron "+var[9].split("/")[0]
                elif var[9]=='':
                    loc = "Exon "+var[8].split("/")[0]

                whitelist = whiteAnnot(p1l, gene, whitelist_file)

                vaf = str(round(100*float(af),2))+"%"

                col = "w3-pale-green w3-hover-green var"

                if float(af)<(LoD/100):
                    col = "w3-pale-red w3-hover-red var"
                    filters = '< LoD'
                if whitelist!="":
                    col = "w3-pale-blue w3-hover-blue var"

                if filters != 'popAF':
                    if filters != 'NM5.25':
                        if int(vd)>3:
                            if keep_csq>0:
                                q=q+1
                                html = html.replace(x,x*2)
                                html = html.replace("##color##",col,1)
                                html = html.replace("##sample_name##",sample_name,1)
                                html = html.replace("##gene##",gene,1)
                                html = html.replace("##dna_pos##",c,1)
                                html = html.replace("##prot_pos1l##",p1l,1)
                                html = html.replace("##prot_pos##",p3l,1)
                                html = html.replace("##localization##",loc,1)
                                html = html.replace("##vaf##",vaf,1)
                                html = html.replace("##csq##",csq,1)
                                html = html.replace("##filters##",filters,1)
                                html = html.replace("##clinvar##",clinvar,1)
                                html = html.replace("##white##",whitelist,1)
                                html = html.replace("##vd##",vd,1)
                                html = html.replace("##dp##",dp,1)
                                html = html.replace("##chrom##",chrom,1)
                                html = html.replace("##pos##",str(pos),1)
                                html = html.replace("##ref##",ref,1)
                                html = html.replace("##alt##",alt,1)
                                html = html.replace("##nm##",nm,1)


    nvar = str(q)
    html = html.replace("##number##",nvar,1)

    html = html.replace(x,"")
    html = html.replace("##sample_name##",sample_name)

    with open(output_report, "w") as FH_out:
        FH_out.write(html)


# create run.html -----------------

run = args.template
run_report = args.output+"/run.html"
file = codecs.open(run, "r", "utf-8").read()

for sample_name in samples:
    file = file.replace('<!--ADD RUN-->','<!--ADD RUN-->'+'\n'+'\t'+'\t'+'\t'+'<!--ADD RUN-->')
    torep = "<tr url="+"\""+sample_name+".html\" "+"""style="cursor: pointer;"><td>"""+sample_name+"</td></tr>"
    file = file.replace('<!--ADD RUN-->',torep,1)

with open(run_report, "w") as FH_out:
    FH_out.write(file)
