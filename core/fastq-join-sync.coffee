###
fastq-join.coffee
(Based on Fastq-join 1.01 by Expression Analysis / Erik Aronesty)

Copyright (c) 2015 Jeongbin Park

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
###

dev = 0
joins = []
#un1s = []
#un2s = []

class _line
    constructor: ->
        this.s = ''
        this.n = 0

class _fq
    constructor: ->
        this.id = new _line
        this.seq = new _line
        this.com = new _line
        this.qual = new _line

comp = (l) ->
    for i in [0...l.length] by 1
        if l[i] == 'A'
            l[i] = 'T'
        else if l[i] == 'T'
            l[i] = 'A'
        else if l[i] == 'G'
            l[i] = 'C'
        else if l[i] == 'C'
            l[i] = 'G'
        else if l[i] == 'a'
            l[i] = 't'
        else if l[i] == 't'
            l[i] = 'a'
        else if l[i] == 'g'
            l[i] = 'c'
        else if l[i] == 'c'
            l[i] = 'g'
    return l

revcomp = (fq) ->
    fq.seq.s = comp(fq.seq.s).reverse() # Non-unicode aware!
    fq.qual.s = fq.qual.s.reverse() # Again, non-unicode aware!
    return fq

min = (a, b) ->
    if a > b
        return b
    else
        return a

max = (a, b) ->
    if a < b
        return b
    else
        return a

hd = (a, s, b, n) ->
    d = 0
    i = 0
    while n > 0
        if a[s+i] != b[i]
            d++
        n--
        i++
    return d+n

# Third mate file is not supported yet
run_fastq_join = (files, donecallback, pgcallback, chunkcallback, joinsonly = false, mino = 6, pctdiff = 8) ->
    dev = 0
    joincnt = 0
    tlen = 0
    tlensq = 0
    nline = [0, 0]
    joins = []
    un1s = []
    un2s = []
    fq = [new _fq, new _fq]
    nrec = 0
    prevfpos = 0
    steprec = 50

    readers = [0, 0]
    for nfin in [0..1] by 1
        entries = files[nfin].name.split('.')
        if entries[entries.length-1] == 'gz'
            gzipped = 1
        readers[nfin] = new jbfilereader(files[nfin], gzipped)

    read_line = (nfin, l, opt) ->
        readers[nfin].readline( (line) -> 
            if readers[0].eof or line == ""
                readers[1].readline( (line2) ->
                    if readers[1].eof or line2 == ""
                        post_loop()
                    else
                        self.postMessage({msgtype: 1, error: "# of rows in mate file doesn't match primary file!"})
                )
            else if readers[1].eof and !readers[0].eof
                self.postMessage({msgtype: 1, error: "# of rows in mate file doesn't match primary file!"})
            else
                l.s = line.split("") # To array, because strings are immutable, note that this is not unicode aware!
                l.n = l.s.length
                opt++
                read_fq(opt)
        )
        return

    post_loop = ->
        chunkcallback(joins)
        pgcallback(100)
        donecallback()
        return
    
    read_fq = (opt) ->
        for nfin in [0..1] by 1
            if opt == 0+nfin*4
                read_line(nfin, fq[nfin].id, opt)
                return
            else if opt == 1+nfin*4
                read_line(nfin, fq[nfin].seq, opt)
                return
            else if opt == 2+nfin*4
                read_line(nfin, fq[nfin].com, opt)
                return
            else if opt == 3+nfin*4
                read_line(nfin, fq[nfin].qual, opt)
                return
        if opt == 8
            #if fq[nfin].id.s[0] == '>'
            #    alert('Input is FASTA!')
            #    return 1
            if (fq[0].id.n == 0 or (fq[0].id.n != 0 and fq[0].id.s[0] != '@')) or
               (fq[1].id.n == 0 or (fq[1].id.n != 0 and fq[1].id.s[0] != '@'))
                self.postMessage({msgtype: 1, error: 'Input file is not FASTQ!'})
                return 1
            main_loop()
            return
        return

    main_loop = ->
        nrec++
        rc = revcomp(fq[1])
        maxo = min(fq[0].seq.n, rc.seq.n)
        bestscore = 2147483647 # INT_MAX of C
        besto = -1
        for i in [mino..maxo] by 1
            mind = (pctdiff * i) // 100
            d = hd(fq[0].seq.s, fq[0].seq.n-i, rc.seq.s, i)
            if d <= mind
                score = (1000*(d*d+1)) // i
                if (score < bestscore)
                    bestscore = score
                    besto = i
        hasex = 0
        olen = besto - hasex
        if besto > 0
            joincnt++

            tlen += olen
            tlensq += olen**2

            #joins.push(fq[0].id.s.join("") + "\n")
            for i in [0...besto] by 1
                li = fq[0].seq.n-besto+i
                ri = i
                if fq[0].seq.s[li] == rc.seq.s[ri]
                    fq[0].qual.s[li] = String.fromCharCode(max(fq[0].qual.s[li].charCodeAt(0), rc.qual.s[ri].charCodeAt(0)))
                else
                    if fq[0].qual.s[li].charCodeAt(0) > rc.qual.s[ri].charCodeAt(0)
                        fq[0].qual.s[li] = String.fromCharCode(33+min(fq[0].qual.s[li].charCodeAt(0),max(fq[0].qual.s[li].charCodeAt(0)-rc.qual.s[ri].charCodeAt(0),3)))
                    else
                        fq[0].seq.s[li] = rc.seq.s[ri]
                        fq[0].qual.s[li] = String.fromCharCode(33+min(rc.qual.s[ri].charCodeAt(0),max(rc.qual.s[ri].charCodeAt(0)-fq[0].qual.s[li].charCodeAt(0),3)))
            joins.push(fq[0].seq.s.join("") + rc.seq.s[besto..].join("") + "\n")
            #joins.push(fq[0].com.s.join("") + "\n")
            #joins.push(fq[0].qual.s.join("") + rc.qual.s[besto..].join("") + "\n")
        #else if !joinsonly
        #    un1s.push(fq[0].id.s.join("") + "\n")
        #    un1s.push(fq[0].seq.s.join("") + "\n")
        #    un1s.push(fq[0].com.s.join("") + "\n")
        #    un1s.push(fq[0].qual.s.join("") + "\n")
        #    un2s.push(fq[1].id.s.join("") + "\n")
        #    un2s.push(fq[1].seq.s.join("") + "\n")
        #    un2s.push(fq[1].com.s.join("") + "\n")
        #    un2s.push(fq[1].qual.s.join("") + "\n")
        if nrec %% steprec == 0
            if prevfpos != readers[0].fpos
                pgcallback(readers[0].fpos*100/files[0].size)
                chunkcallback(joins)
                joins = []
                prevfpos = readers[0].fpos
            setTimeout( ( -> read_fq(0)), 0)
        else
            read_fq(0)
        return
    setTimeout( ( -> read_fq(0)), 0)
    return
