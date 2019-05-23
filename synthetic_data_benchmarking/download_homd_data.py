import os

lines=[]
os.system('wget http://www.homd.org/ftp/HOMD_annotated_genomes/')
with open('index.html') as f:
    for line in f:
        try:
            lines.append(line.split('href="')[1].split('/">')[0])
        except:
            continue

lines=lines[2:-1]
os.system('rm index.html')

for line in lines:
  os.system('wget "http://www.homd.org/ftp/HOMD_annotated_genomes/%s/"'%line)
  lines2=[]
  with open('index.html') as f:
    for line2 in f:
      if 'href' in line2:
        lines2.append(line2.split('href')[1])
  lines2=[x.split('</a>')[0].split('>')[1] for i,x in enumerate(lines2) if i<6]
  fullSeq=lines2[2]
  genes=lines2[-1]
  os.system('axel -a -o %s http://www.homd.org/ftp/HOMD_annotated_genomes/%s/%s'%(line+'_fullSeq',line,fullSeq))
  os.system('axel -a http://www.homd.org/ftp/HOMD_annotated_genomes/%s/%s -o %s'%(line,genes,line+'_genes'))
  os.system('rm index.html')
