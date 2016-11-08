###
jbfilereader.coffee

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

class _inflater
    constructor: ->
        window_bits = 15
        enable_windows_gzip = 32
        @inflater_pako = new pako.Inflate({to: 'string', chunkSize: 16384, windowBits: window_bits|enable_windows_gzip})
        @inflating_buffer = ''
        that = this
        @inflater_pako.onData = (chunk) ->
            that.inflating_buffer += chunk
            return
        return

    decompress: (chunk, islastchunk) ->
        @inflating_buffer = ''
        @inflater_pako.push(chunk, islastchunk)
        @ended = @inflater_pako.ended
        @strm = @inflater_pako.strm
        return @inflating_buffer

class jbfilereadersync
    constructor: (@file, @gzipped) ->
        @buffer = ''
        @filesize = @file.size
        @chunksize = 1024 * 512
        @reader = new FileReaderSync()

        if @gzipped
            @inflater = new _inflater()

        @islastchunk = false

        @fpos = 0
        @endpos = 0

        @eof = false

        return

    _getchunk: ->
        if @fpos + @chunksize >= @filesize
            @endpos = @filesize
            @islastchunk = true
        else
            @endpos = @fpos + @chunksize
        blob = @file.slice(@fpos, @endpos)

        @fpos += @endpos - @fpos
        if @gzipped
            raw_array = new Uint8Array(@reader.readAsArrayBuffer(blob))
            s = @inflater.decompress(raw_array, @islastchunk)
            if s
                if @inflater.ended and (@inflater.strm.avail_in or !@islastchunk)
                    # Non-standard gzip. See http://www.gzip.org/#faq8
                    remaining_bytes = @inflater.strm.avail_in
                    rel_pos = 0
                    while raw_array[raw_array.byteLength-remaining_bytes+rel_pos] == 0
                        rel_pos++
                    @fpos -= remaining_bytes-rel_pos
                    @inflater = new _inflater() # Renew Pako
            else
                throw 'Something wrong with the gzipped file!'
        else
            s = @reader.readAsText(blob)
        return s

    readline: () ->
        if @eof
            return ""

        lfpos = @buffer.indexOf("\n")
        while lfpos == -1
            if @fpos >= @filesize
                result = @buffer
                @buffer = ""
                @fpos = @filesize
                @eof = true
                return result
            @buffer += @_getchunk()
            lfpos = @buffer.indexOf("\n")

        if @buffer[lfpos-1] == "\r"
            result = @buffer[...lfpos-1]
        else
            result = @buffer[...lfpos]
        @buffer = @buffer[lfpos+1...]
        return result

class jbfilereader
    constructor: (@file, @gzipped) ->
        @buffer = ''
        @filesize = @file.size
        @chunksize = 1024 * 512
        @reader = new FileReader()

        if @gzipped
            @inflater = new _inflater()

        @islastchunk = false

        @fpos = 0
        @endpos = 0

        @eof = false

        return

    _readblob: (blob) ->
        that = this
        readpromise = new Promise( (resolve, reject) ->
            that.reader.onload = (e) ->
                resolve(e.target.result)
            that.reader.onerror = ->
                reject()
            return
        )
        if @gzipped
            @reader.readAsArrayBuffer(blob)
        else
            @reader.readAsText(blob)
        return readpromise

    _getchunk: ->
        if @fpos + @chunksize >= @filesize
            @endpos = @filesize
            @islastchunk = true
        else
            @endpos = @fpos + @chunksize
        blob = @file.slice(@fpos, @endpos)

        that = this
        chunkpromise = new Promise( (resolve, reject) ->
            readpromise = that._readblob(blob)
            readpromise.then( (s) ->
                that.fpos += that.endpos - that.fpos
                if that.gzipped
                    raw_array = new Uint8Array(s)
                    s = that.inflater.decompress(raw_array, that.islastchunk)
                    if s
                        if that.inflater.ended and (that.inflater.strm.avail_in or !that.islastchunk)
                            # Non-standard gzip. See http://www.gzip.org/#faq8
                            remaining_bytes = that.inflater.strm.avail_in
                            rel_pos = 0
                            while raw_array[raw_array.byteLength-remaining_bytes+rel_pos] == 0
                                rel_pos++
                            that.fpos -= remaining_bytes-rel_pos
                            that.inflater = new _inflater() # Renew Pako
                    else
                        alert('Something wrong with the gzipped file!')
                        reject()
                that.buffer += s
                resolve()
                return
            ).catch( (s) ->
                alert('Something wrong while reading file!')
                reject()
                return
            )
        )
        return chunkpromise

    readline: (callback) ->
        if @eof
            callback("")

        lfpos = @buffer.indexOf("\n")
        if lfpos == -1
            if @fpos >= @filesize
                result = @buffer
                @buffer = ""
                @eof = true
                callback(result)
            else
                that = this
                datapromise = @_getchunk()
                datapromise.then( -> that.readline(callback) )
        else
            if @buffer[lfpos] == "\r"
                result = @buffer[...lfpos-1]
            else
                result = @buffer[...lfpos]
            @buffer = @buffer[lfpos+1...]
            callback(result)
        return
