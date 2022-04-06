import difflib
import sys

diff = difflib.HtmlDiff().make_file("AAGGCCTT","GAGGCCGG")

f=open('test.html', 'w')
f.write(diff)
f.close()