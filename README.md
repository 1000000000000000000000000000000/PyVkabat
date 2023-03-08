# PyVkabat
PyVkabat is a bioinformatics/computational chemistry tool used to calculated secondary structure variability at each residue position in a protein of interest. To reduce run time, each of the requests are made in parallel.

<br>
Vkabat = k * N / n1

where:
<ul><li>k is the number of secondary structure classes predicted for a residue</li>
    <li>N is the total number of predictions made for a given residue</li>
    <li>n1 is how often the most predicted secondary structure is observed for a residue</li></ul><br>
    
PyVkabat submits the user's protein sequence of interest to 15 different algorithms:
<ul>
  <li>gor1</li>
  <li>gor3</li>
  <li>dpm</li>
  <li>sopm</li>
  <li>predator</li>
  <li>hnn</li>
  <li>dsc</li>
  <li>mlrc</li>
  <li>yaspin</li>
  <li>jpred</li>
  <li>phd</li>
  <li>prof</li>
  <li>sspro</li>
  <li>jnet</li>
  <li>psipred</li>
</ul>
