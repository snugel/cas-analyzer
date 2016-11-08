/*
 * needle.js
 *
 * Author: Jeongbin Park (pjb7687 at gmail.com)
 * Based on source code of 'needle' in EMBOSS software suite (http://emboss.sourceforge.net/).
 *
 * GNU General Public License
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

    function needle(a, b, gapopen, gapextend, endgapopen, endgapextend) {
            function E_FPEQ(a, b, e) { return (((b - e) < a) && (a < (b + e))); }

            var sub = [[5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -4], [-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5], [-4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, -4], [-4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -4], [-4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -4], [1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, 1], [1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -4], [-4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, 1], [-4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, 1], [1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -4], [-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -1], [-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -4], [-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -1], [-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -1], [-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2], [-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5]]; // EDNAFULL
            var base_to_idx = {'A': 0, 'C': 3, 'B': 10, 'D': 13, 'G': 2, 'H': 12, 'K': 8, 'M': 9, 'N': 14, 'S': 4, 'R': 6, 'U': 15, 'T': 1, 'W': 5, 'V': 11, 'Y': 7}; // for EDNAFULL
        var idx_to_base = {0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: 'S', 5: 'W', 6: 'R', 7: 'Y', 8: 'K', 9: 'M', 10: 'B', 11: 'V', 12: 'H', 13: 'D', 14: 'N', 15: 'U'}; // for EDNAFULL

        var ix = [], iy = [], m = [];
        var ixp, iyp, mp;

        var compass = [];

        var i, j;
        var cursor, cursorp;

        var p = []; q = [];

        for(i=0; i<a.length; i++) p.push(base_to_idx[a[i]]);
        for(i=0; i<b.length; i++) q.push(base_to_idx[b[i]]);

        var lena = a.length;
        var lenb = b.length;

        var match = sub[p[0]][q[0]];
        var score;

        var aln_a = [], aln_b = [], aln_r = [];

        ix[0] = -endgapopen-gapopen;
        iy[0] = -endgapopen-gapopen;
        m[0] = match;

        for (ypos=1; ypos<lena; ypos++) {
            match = sub[p[ypos]][q[0]];
            cursor = ypos * lenb;
            cursorp = (ypos - 1) * lenb;

            testog = m[cursorp] - gapopen;
            testeg = iy[cursorp] - gapextend;

            if (testog >= testeg)
                iy[cursor] = testog;
            else
                iy[cursor] = testeg;
            m[cursor] = match - (endgapopen + (ypos - 1) * endgapextend);
            ix[cursor] = -endgapopen - ypos * endgapextend - gapopen;
        }
        ix[cursor] -= endgapopen;
        ix[cursor] += gapopen;

        cursor = 0;

        for (xpos=1; xpos<lenb; xpos++) {
            match = sub[p[0]][q[xpos]];
            cursor = xpos;
            cursorp = xpos - 1;
         
            testog = m[cursorp] - gapopen;
            testeg = ix[cursorp] - gapextend;
            
            if(testog >= testeg)
                ix[cursor] = testog;
            else
                ix[cursor] = testeg;

            m[cursor] = match - (endgapopen + (xpos - 1) * endgapextend);
            iy[cursor] = -endgapopen - xpos * endgapextend - gapopen;
        }
        iy[cursor] -= endgapopen;
        iy[cursor] += gapopen;

        xpos = 1;
        while (xpos != lenb) {
            ypos = 1;
            bconvcode = q[xpos];

            cursorp = xpos-1;
            cursor = xpos++;

            while (ypos < lena) {
                match = sub[p[ypos++]][bconvcode];
                cursor += lenb;
                
                mp = m[cursorp];
                ixp = ix[cursorp];
                iyp = iy[cursorp];

                if (mp > ixp && mp > iyp)
                    m[cursor] = mp + match;
                else if (ixp > iyp)
                    m[cursor] = ixp + match;
                else
                    m[cursor] = iyp + match;

                cursorp += 1;
                if (xpos == lenb) {
                    testog = m[cursorp] - endgapopen;
                    testeg = iy[cursorp] - endgapextend;
                } else {
                    testog = m[cursorp];
                    if (testog < ix[cursorp])
                        testog = ix[cursorp];

                    testog -= gapopen;
                    testeg = iy[cursorp] - gapextend;
                }
                if (testog > testeg)
                    iy[cursor] = testog;
                else
                    iy[cursor] = testeg;

                cursorp += lenb;
                
                cursorp -= 1;
                if (ypos == lena) {
                    testog = m[cursorp] - endgapopen;
                    testeg= ix[cursorp] - endgapextend;
                } else {
                    testog = m[cursorp];
                    if (testog < iy[cursorp])
                        testog = iy[cursorp];
                    testog -= gapopen;
                    testeg = ix[cursorp] - gapextend;
                }
                if (testog > testeg)
                    ix[cursor] = testog;
                else
                    ix[cursor] = testeg;
            }
        }

        score = -32767; // INT_MIN
        start1 = lena - 1;
        start2 = lenb - 1;

        cursor = lena * lenb - 1;
        if (m[cursor] > ix[cursor] && m[cursor] > iy[cursor])
            score = m[cursor];
        else if (ix[cursor] > iy[cursor])
            score = ix[cursor];
        else
            score = iy[cursor];

        cursorp = 0;
        cursor = 1;
        
        eps = 1.192e-6;

        ypos = start1;
        xpos = start2;

        while (xpos>=0 && ypos>=0) {
            cursor = ypos*lenb+xpos;
            mp = m[cursor];

            if(cursorp == 1 && E_FPEQ((ypos==0||(ypos==lena-1)?endgapextend:gapextend), (ix[cursor]-ix[cursor+1]), eps))
            {
                compass[cursor] = 1;
                xpos--;
            }
            else if(cursorp== 2 && E_FPEQ((xpos==0||(xpos==lenb-1)?endgapextend:gapextend), (iy[cursor]-iy[cursor+lenb]), eps))
            {
                compass[cursor] = 2;
                ypos--;
            }
            else if(mp >= ix[cursor] && mp>= iy[cursor])
            {

                if(cursorp == 1 && E_FPEQ(mp,ix[cursor],eps))
                {
                    compass[cursor] = 1;
                    xpos--;
                }
                else if(cursorp == 2 && E_FPEQ(mp,iy[cursor],eps))
                {
                    compass[cursor] = 2;
                    ypos--;
                }
                else
                {
                    compass[cursor] = 0;
                    ypos--;
                    xpos--;
                }

            }
            else if(ix[cursor]>=iy[cursor] && xpos>-1)
            {
                compass[cursor] = 1;
                xpos--;
            }
            else if(ypos>-1)
            {
                compass[cursor] = 2;
                ypos--;
            }
            else
            {
                alert("Needle: Something is seriously wrong in the traceback algorithm");
                return -1;
            }
            cursorp = compass[cursor];
        }

        for (i=lenb-1; i>start2;) {
            aln_b.push(idx_to_base[q[i--]]);
            aln_a.push('-');
            aln_r.push(' ');
        }
        for (j=lena-1; j>start1;) {
            aln_a.push(idx_to_base[p[j--]]);
            aln_b.push('-');
            aln_r.push(' ');
        }

        while (start2 >= 0 && start1 >= 0)
        {
            cursor = start1 * lenb + start2;
            if(!compass[cursor]) /* diagonal */
            {
                b1 = p[start1--]; b2 = q[start2--];
                aln_a.push(idx_to_base[b1]);
                aln_b.push(idx_to_base[b2]);
                if (b1 == b2)
                    aln_r.push('|');
                else
                    aln_r.push('.');
                continue;
            }
            else if(compass[cursor] == 1) /* Left, gap(s) in vertical */
            {
                aln_a.push('-');
                aln_b.push(idx_to_base[q[start2--]]);
                aln_r.push(' ');
                continue;
            }
            else if(compass[cursor] == 2) /* Down, gap(s) in horizontal */
            {
                aln_a.push(idx_to_base[p[start1--]]);
                aln_b.push('-');
                aln_r.push(' ');
                continue;
            } else {
                alert("Needle: Walk Error in NW");
                return -1;
            }
        }

        for (;start2>=0;start2--)
        {
            aln_b.push(idx_to_base[q[start2]]);
            aln_a.push('-');
            aln_r.push(' ');
        }

        for (;start1>=0;start1--)
        {
            aln_a.push(idx_to_base[p[start1]]);
            aln_b.push('-');
            aln_r.push(' ');
        }

        aln_a = aln_a.reverse().join("");
        aln_b = aln_b.reverse().join("");
        aln_r = aln_r.reverse().join(""); //.replace(/ /gi, "&nbsp;");

        return [aln_a, aln_r, aln_b];
    }
