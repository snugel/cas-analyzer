importScripts('core/pako_inflate.min.js');
importScripts('core/jbfilereader.js');
importScripts('core/fastq-join.js');
importScripts('core/needle.js');
importScripts('core/miseq-analyzer.js');

var rgen_type, seq_RGEN, seq_RGEN2, seq_wt, seq_hdr, end_range, filt_n, filt_r, files;

onmessage = function(msg) {
	if (msg.data.msgtype == 0) {
		rgen_type = msg.data.rgen_type;
		seq_RGEN = msg.data.seq_RGEN;
		seq_RGEN2 = msg.data.seq_RGEN2;
		seq_wt = msg.data.seq_wt;
		seq_hdr = msg.data.seq_hdr;
		end_range = msg.data.end_range;
		filt_n = msg.data.filt_n;
		filt_r = msg.data.filt_r;
		files = msg.data.files;
    
		if (rgen_type < 2) {
			pattern = new RegExp("(" + seq_RGEN + ")|(" + revcompstr(seq_RGEN) + ")");
			m = pattern.exec(seq_wt);
			if (m) {
				if (rgen_type === 0) {
					if (m[1]) {
						bp = m.index + seq_RGEN.length - 3;
					} else {
						bp = m.index + 3;
					}
				} else if (rgen_type === 1) {
					if (m[1]) {
						bp = m.index + 21;
					} else {
						bp = m.index + seq_RGEN.length - 16;
					}
				}
			} else {
				self.postMessage({msgtype: 1, error: "Couldn't find RGEN site in WT sequence!"});
				return;
			}
		} else {
			if (rgen_type === 2) {
				pattern_1 = new RegExp("(" + seq_RGEN + ")");
				pattern_2 = new RegExp("(" + seq_RGEN2 + ")");
			} else {
				pattern_1 = new RegExp("(" + seq_RGEN + ")|(" + revcompstr(seq_RGEN) + ")");
				pattern_2 = new RegExp("(" + seq_RGEN2 + ")|(" + revcompstr(seq_RGEN2) + ")");
			}
			m_1 = pattern_1.exec(seq_wt);
			m_2 = pattern_2.exec(seq_wt);
			if (m_1 && m_2) {
				if (m_1.index > m_2.index) {
					self.postMessage({msgtype: 1, error: "Position of left site is bigger than that of right site!"});
					return;
				}
				if (rgen_type < 4) {
					bp = Math.round((m_1.index + m_2.index + seq_RGEN2.length) / 2);
				} else if (rgen_type === 4) {
					if (m_1[1]) {
						bp_1 = m.index + seq_RGEN.length - 3;
					} else {
						bp_1 = m.index + 3;
					}
					if (m_2[1]) {
						bp_2 = m.index + seq_RGEN2.length - 3;
					} else {
						bp_2 = m.index + 3;
					}
					bp = Math.round((bp_1 + bp_2) / 2);
				} else {
					if (m_1[1]) {
						bp_1 = m.index + 21;
					} else {
						bp_1 = m.index + seq_RGEN.length - 16;
					}
					if (m_2[1]) {
						bp_2 = m.index + 21;
					} else {
						bp_2 = m.index + seq_RGEN2.length - 16;
					}
					bp = Math.round((bp_1 + bp_2) / 2);
				}
			} else {
				self.postMessage({msgtype: 1, error: "Couldn't find both sites in WT sequence!"});
				return;
			}
		}
		start_pos = bp - end_range;
		end_pos = bp + end_range;
		if (start_pos < 0) {
			start_pos = 0;
		}
		if (end_pos > seq_wt.length) {
			end_pos = seq_wt.length;
		}
		s_seq = seq_wt.slice(bp - filt_r, bp + filt_r);
		seq_range = seq_wt.slice(start_pos, end_pos);
		pri_for = seq_range.slice(0, pri_len);
		pri_back = seq_range.slice(-pri_len);
		pri_for_patterns = [];
		pri_back_patterns = [];
		for (i = 0; i < pri_len; i++) {
			pri_for_patterns.push(pri_for.slice(0, i) + '[AGCT]' + pri_for.slice(i + 1));
			pri_back_patterns.push(pri_back.slice(0, i) + '[AGCT]' + pri_back.slice(i + 1));
		}
		min = function(a, b) {
			if (a < b) {
				return a;
			} else {
				return b;
			}
		};
		max = function(a, b) {
			if (a > b) {
				return a;
			} else {
				return b;
			}
		};
		seq_fancy_wt = seq_wt.slice(0, start_pos);
		seq_fancy_wt += '<span class="primerseq">' + seq_wt.slice(start_pos, min(start_pos + pri_len, bp - filt_r)) + '</span>';
		seq_fancy_wt += seq_wt.slice(start_pos + pri_len, min(m.index, bp - filt_r));
		seq_fancy_wt += '<span class="rgenseq">' + seq_wt.slice(m.index, bp - filt_r) + '</span>';
		seq_fancy_wt += '<span class="coreseq">' + seq_wt.slice(bp - filt_r, bp + filt_r) + '</span>';
		seq_fancy_wt += '<span class="rgenseq">' + seq_wt.slice(bp + filt_r, m.index + seq_RGEN.length) + '</span>';
		seq_fancy_wt += seq_wt.slice(max(m.index + seq_RGEN.length, bp + filt_r), end_pos - pri_len);
		seq_fancy_wt += '<span class="primerseq">' + seq_wt.slice(max(end_pos - pri_len, bp + filt_r), end_pos) + '</span>';
		seq_fancy_wt += seq_wt.slice(end_pos);

		set_primer_patterns(pri_for_patterns, pri_back_patterns);

		self.postMessage({msgtype: 0, seq_fancy_wt: seq_fancy_wt});    
	} else if (msg.data.msgtype == 1 || msg.data.msgtype == 2) {
		empty_cache();
		pgcallback1 = function(p) { self.postMessage({msgtype: 2, progress: p}); };
		pgcallback2 = function(p) { self.postMessage({msgtype: 3, progress: p}); };
		if (msg.data.msgtype == 1)
			run_fastq_join(msg.data.files, pgcallback1, process_chunk, true);
		else
			parse_file(msg.data.file, pgcallback1, process_chunk);
		data = run_cas_analyser(seq_range, seq_hdr, filt_n, filt_r, pgcallback2);
	        self.postMessage({msgtype: 4, results: data});
	}
}
