Summary of project Statement of Work document

First test of pandoc for tex -> docx


After installing pandoc, navigate to directory, and type "chcp 65001" to ensure  the encoding is UTF-8.
Then type:  
pandoc SomeFile.tex -s -o SomeFile.docx  

-s flag means create standalone file   
-o SomeFile.*** flag means output to specific file  
Learned from:  
https://pandoc.org/getting-started.html  

Doesn't seem to catch/handle headers from tex->docx  