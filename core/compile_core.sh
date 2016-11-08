cp jbfilereader/jbfilereader.coffee .
cp pako/dist/pako_inflate.min.js .
cp needle_js/needle.js .
coffee --bare --compile jbfilereader.coffee
coffee --bare --compile fastq-join.coffee
coffee --bare --compile miseq-analyzer.coffee
