###
interface.coffee

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

plots = []
files = []
joinedfile = 0
singleendfile = 0

seq_wt = ''
seq_fancy_wt = ''
seq_RGEN = ''
seq_RGEN2 = ''
filt_r = 0
filt_n = 0
end_range = 0
joins_length = 0
joins = []

skip_progress = 0

seq_hdr = ''
rgen_type = 0 # 0 for Cas9, 1 for Cpf1, 2 for TALEN/ZFN/dCas9-FokI, 3 for paired Cas9, 4 for paired Cpf1
cate = -1
_results = []

worker = undefined

#prevsinglefile = ''
#prevfile1 = ''
#prevfile2 = ''

update_progress_fastqjoin = (progress) ->
    $('#progress_fastqjoin').css('width', progress + '%')

$('#submit').click ->
    if typeof(Worker) != "undefined"
        worker = new Worker("/static/cas-analyzer/worker.js")
        window.cas_analyzer_worker = worker # Keep reference of worker so that it won't garbage collected
    else
        alert("Your browser currently does not support features used by Cas-Analyzer. Please consider upgrading.")
    seq_wt = new String($("#fullseq").val()).toUpperCase().replace(/\s/g,'')
    seq_hdr = new String($("#hdrseq").val()).toUpperCase().replace(/\s/g,'')

    if seq_hdr
        n = seq_wt.search(seq_hdr[...10])+10
        lcnt = 10
        while seq_wt[n+lcnt-10] == seq_hdr[lcnt]
            lcnt += 1

        n = seq_wt.search(seq_hdr[-10...])-1
        rcnt = 11
        while seq_wt[n-rcnt+11] == seq_hdr[-rcnt...][0]
            rcnt += 1
    
        while seq_hdr.length-lcnt-rcnt < 10
            lcnt -= 1
            rcnt -= 1

        seq_hdr = seq_hdr[lcnt...seq_hdr.length+1-rcnt]

    # rgen_type: 0 for Cas9, 1 for Cpf1, 2 for TALEN/ZFN, 3 for dCas9-FokI, 4 for paired Cas9, 5 for paired Cpf1
    seq_RGEN = new String($("#rgenseq").val()).toUpperCase().replace(/\s/g,'')
    if $("#nuctype").val() == 0 and $("#nucleases").val() == 6
        rgen_type = 1
    else if $("#nuctype").val() == 1
        seq_RGEN2 = new String($("#rgenseq-right").val()).toUpperCase().replace(/\s/g,'')
        if $("#nucleases").val() == 9
            if $("#pams").val() == 6
                rgen_type = 5
            else
                rgen_type = 4
        else if $("#nucleases").val() == 10
            rgen_type = 3
        else
            rgen_type = 2

    filt_n = parseInt($("#nval").val())
    filt_r = parseInt($("#rval").val())
    if $('#chklr').is(":checked")
        end_range = seq_wt.length
    else
        end_range = parseInt($("#lrval").val())

    if isNaN(end_range)
        end_range = 70
    if isNaN(filt_n)
        filt_n = 1
    if $('#chkr').is(":checked")
        if isNaN(filt_r)
            filt_r = 5
    else
        filt_r = 0

    fileopt = parseInt($('#optfile').val())
    if fileopt == 0
        files[0] = document.getElementById("file1").files[0]
        files[1] = document.getElementById("file2").files[0]
        if !files[0]
            alert('Please set file R1!')
            return
        if !files[1]
            alert('Please set file R2!')
            return
    else if fileopt == 1
        singleendfile = document.getElementById("file3").files[0]
        if !singleendfile
            alert('Please set single read file!')
            return
    if !seq_wt
        alert('Please input wildtype sequence!')
        return
    if !seq_RGEN
        alert('Please input RGEN sequence!')
        return
    if end_range <= 12
        alert('R should be at least 12!')
        return
    if fileopt == 0
        $('.file-loading-results-header').html('<th colspan="3">File1</th><th colspan="3">File2</th>')
        $('.joinedfn1').attr('colspan', '3')
        $('.joinedfn2').show()
        $('.joinedfn1').html(files[0].name)
        $('.joinedfn2').html(files[1].name)
    else
        $('.file-loading-results-header').html('<th colspan="6">File name</th>')
        $('.joinedfn1').attr('colspan', '6')
        $('.joinedfn2').hide()
        $('.singlefn').html(singleendfile.name)
    $('.input-rgenseq').html(seq_RGEN)
    $('.input-rval').html(String(filt_r))
    $('.input-nval').html(String(filt_n))
    if $('#chklr').is(':checked')
        $('.input-lrval').html('Using both ends')
    else
        $('.input-lrval').html(String(end_range))
    data = {}
    data.msgtype = 0
    data.seq_wt = seq_wt
    data.rgen_type = rgen_type
    data.seq_RGEN = seq_RGEN
    data.seq_RGEN2 = seq_RGEN2
    data.seq_hdr = seq_hdr
    data.end_range = end_range
    data.filt_n = filt_n
    data.filt_r = filt_r
    data.files = files

    worker.onmessage = (msg) ->
        if msg.data.msgtype == 0
            $('.input-fullseq').html(msg.data.seq_fancy_wt)
            location.hash = "!progress"
            if fileopt == 0
                $("#upper_progress").text("Fastq-join")
                setTimeout( ->
                    #if prevsinglefile == '' and prevfile1 == (files[0].name+files[0].size) and prevfile2 == (files[1].name+files[1].size)
                    #    worker.postMessage({msgtype: 3})
                    #else
                    #    prevfile1 = files[0].name+files[0].size
                    #    prevfile2 = files[1].name+files[1].size
                    #    prevsinglefile = ''
                    #    worker.postMessage({msgtype: 1, files: files})
                    #return
                    worker.postMessage({msgtype: 1, files: files})
                , 0)
            else if fileopt == 1
                $("#upper_progress").text("File loading")
                setTimeout( ->
                    #if prevsinglefile == singleendfile.name+singleendfile.size and prevfile1 == '' and prevfile2 == ''
                    #    worker.postMessage({msgtype: 3})
                    #else
                    #    prevsinglefile = singleendfile.name+singleendfile.size
                    #    prevfile1 = ''
                    #    prevfile2 = ''
                    #    worker.postMessage({msgtype: 2, file: singleendfile})
                    #return
                    worker.postMessage({msgtype: 2, file: singleendfile})
                , 0)
        else if msg.data.msgtype == 1
            alert(msg.data.error) # Error
            location.hash = ""
        else if msg.data.msgtype == 2
            setTimeout( ->
                $('#progress_fastqjoin').css('width', msg.data.progress + '%')
                return
            , 0)
        else if msg.data.msgtype == 3
            setTimeout( ->
                $('#progress_analysis').css('width', msg.data.progress + '%')
                return
            , 0)
        else
            all_done(msg.data.results)
    worker.postMessage(data)
    return

$('#btnselect10').click ->
    $('.chkmsas').attr('checked', false)
    for i in [0..9] by 1
        $('#chkmsa'+i).attr('checked', true)
    return

$('#chklr').click ->
    if $(this).is(":checked")
        $('#lrval').attr('disabled', true)
    else
        $('#lrval').attr('disabled', false)
    return

$('#chkr').click ->
    if $(this).is(":checked")
        $('#rval').attr('disabled', false)
    else
        $('#rval').attr('disabled', true)
    return

$('#btnshowall').click( (e) ->
    e.preventDefault()
    $('.seqfilter').removeClass("active")
    $(this).parent().addClass("active")
    cate = -1
    list_data()
)

$('#btnshowws').click( (e) ->
    e.preventDefault()
    $('.seqfilter').removeClass("active")
    $(this).parent().addClass("active")
    cate = 0
    list_data()
)

$('#btnshowins').click( (e) ->
    e.preventDefault()
    $('.seqfilter').removeClass("active")
    $(this).parent().addClass("active")
    cate = 1
    list_data()
)

$('#btnshowdel').click( (e) ->
    e.preventDefault()
    $('.seqfilter').removeClass("active")
    $(this).parent().addClass("active")
    cate = 2
    list_data()
)

$('#showhdr').click( (e) ->
    list_data()
)

list_data = ->
    ishdr = $('#showhdr').is(":checked")
    view_data = []
    for i in [0..._results.length] by 1
        if (cate == -1 or _results[i][6] == cate) and (!ishdr or (ishdr and _results[i][7] > 0))
            view_data.push(_results[i])
    load_page(0, view_data, 10)

all_done = (data) ->
    joins_count = data.joins_length
    totlr_count = data.totlr_count
    totn_count = data.tot_count
    ins_count = data.cnt_ins
    del_count = data.cnt_del
    results = data.table
    il = data.il
    dl = data.dl
    isi = data.is
    dsi = data.ds
    hdr = data.hdr

    $('.result-joined').html(String(joins_count))
    $('.result-totallr').html(String(totlr_count))
    $('.result-totaln').html(String(totn_count))
    $('.result-insertions').html(String(ins_count))
    $('.result-deletions').html(String(del_count))
    $('.result-indelratio').html(String(ins_count+del_count)+" ("+String(((ins_count+del_count)*100.0/totn_count).toFixed(1)) + "%)")
    $('.result-hdrratio').html(String(hdr)+" ("+String((hdr*100.0/totn_count).toFixed(1)) + "%)")
    rows = ''
    if results.length == 0
        rows = '<tr><td colspan=5>No entries found.</td></tr>'
    else
        ins_list = []
        del_list = []
        ws_list = []

    _results = results

    $('#btndownload').off()
    $('#btndownload').on('click', ->
        s = '#ID\tWT Sequence\tRGEN Treated Sequence\tLength\tCount\tType\tHDR\n'
        for i in [0...results.length] by 1
            entry = results[i]
            s += entry[0] + '\t' + entry[1] + '\t' + entry[2] + '\t' + entry[4] + '\t' + entry[5] + '\t'
            if entry[6] == 0
                s += 'WT or Sub'
            else if entry[6] == 1
                s += 'Ins'
            else
                s += 'del'
            s += '\t'
            if entry[7] < -1
                s += 'N/A'
            else if entry[7] == -1
                s += 'X'
            else
                s += 'O'
            s += '\n'
        saveTextData(s, 'result.txt')
        return
    )

    draw_plots = ->
        if ($("#placeholder_is").width())
            load_page(0, results, 10)
            plots.push(
                $.plot("#placeholder_is", [
                    { data: isi, color: "#DE000F" }
                ], {
                    series: {
                        bars: { show: true },
                    },
                    grid: {
                        hoverable: true
                    },
                    xaxis: {
                        #transform: (v) -> return if v == 0 then v else Math.log(v),
                        axisLabel: "Size",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 10
                    },
                    yaxis: {
                        axisLabel: "Count",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 8
                    }
                })
            )
            plots.push(
                $.plot("#placeholder_ds", [
                    { data: dsi, color: "#DE000F" }
                ], {
                    series: {
                        bars: { show: true },
                    },
                    grid: {
                        hoverable: true
                    },
                    xaxis: {
                        #transform: (v) -> return if v == 0 then v else Math.log(v),
                        axisLabel: "Size",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 10
                    },
                    yaxis: {
                        axisLabel: "Count",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 8
                    }
                })
            )
            plots.push(
                $.plot("#placeholder_il", [
                    { data: il, color: "#5482FF" }
                ], {
                    series: {
                        bars: { show: true },
                    },
                    grid: {
                        hoverable: true
                    },
                    xaxis: {
                        axisLabel: "Position",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 10
                    },
                    yaxis: {
                        axisLabel: "Ratio (%)",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 8
                    }
                })
            )
            plots.push(
                $.plot("#placeholder_dl", [
                    { data: dl, color: "#5482FF" }
                ], {
                    series: {
                        bars: { show: true },
                    },
                    grid: {
                        hoverable: true
                    },
                    xaxis: {
                        axisLabel: "Position",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 10
                    },
                    yaxis: {
                        axisLabel: "Ratio (%)",
                        axisLabelUseCanvas: true,
                        axisLabelFontFamily: 'Arial',
                        axisLabelPadding: 8
                    }
                })
            )
            $("<div id='tooltip_plot'></div>").css({
                position: "absolute",
                display: "none",
                border: "1px solid #fdd",
                padding: "2px",
                "background-color": "#fee",
                opacity: 0.80
            }).appendTo("body")
            $(".plotarea").each( ->
                $(this).bind("plothover", (evt, pos, item) ->
                    if (item)
                        x = item.datapoint[0]
                        if this.id == "placeholder_ds" or this.id == "placeholder_is"
                            y = item.datapoint[1]
                        else
                            y = item.datapoint[1].toFixed(3)
                        $("#tooltip_plot").html("(" + x + ", " + y + ")").css({top: item.pageY+5, left: item.pageX+5}).fadeIn(200)
                    else
                        $("#tooltip_plot").hide()
                    return
                )
                return
            )
        else
            setTimeout(draw_plots, 200)
        return
    setTimeout(draw_plots, 200)

    location.hash = "!result"
    skip_progress = 1
    return
clear_plots = ->
    for plot in plots
        plot.shutdown()
    plots = []
window.onload = ->
    location.hash = "!"
    return

window.onhashchange = ->
    $('#progress_fastqjoin').css('width', '0')
    $('#progress_analysis').css('width', '0')
    if location.hash == "#!progress"
        if skip_progress
            skip_progress = 0
            location.hash = "!"
        else
            $("#inputform").hide()
            $("#progress").show()
    else if location.hash == "#!result"
        $("#progress").hide()
        $("#result").show()
    else
        worker.terminate()
        clear_plots()
        $("#progress").hide()
        $("#result").hide()
        $("#inputform").show()
        $('.seqfilter').removeClass("active")
        $('#btnshowall').parent().addClass("active")
    return
