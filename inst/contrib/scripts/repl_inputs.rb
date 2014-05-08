#! /usr/bin/ruby
#
# This is a script to replace all \input{} statements in the Rd documentation
# and copyright statements in file headers. This to prevent duplication of text
# all over the place. In the future we may use the \Sexpr{} macro - from R 2.10.0.
#
# In principle a line like
#
#   % \input{"inst/docs/mqm/limitations.txt"}
#   REPLACE THIS
#   }
# 
# will inject the contents of that file - replacing the text until the next
# closing curly brace. There are some other replacement commands that
# merely replace LaTeX like macros on the next line
#
#   % \mqmcopyright   # default copyright
#
# Usage:
#
#   ruby ./contrib/scripts/repl_inputs.rb inputfile(s)
#
# Example
#
#   ruby ./contrib/scripts/repl_inputs.rb man/*.Rd
#

REPL = 
{ 
  # Keywords for C/C++ headers
  'mqmcopyright1' =>
    " * Copyright (c) 1998-2010, Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman",
  'license' => 
    ' * Published under the terms of the GNU Lesser General Public License 3 (GPL3)',
  # Keywords for R man pages
  'mqmauthors' =>
    "Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \\email{kbroman@biostat.wisc.edu}",
  'dannyauthor' =>
    "Danny Arends \\email{danny.arends@gmail.com}",
  'dannyrutgerauthors' =>
    "Danny Arends \\email{danny.arends@gmail.com} ; Rutger Brouwer",
  'crossobject' =>
    'An object of class \code{cross}. See \code{\link{read.cross}} for details.',
  'mqmscanobject' =>
    'An object returned by \code{mqmscan}, including cofactors and QTL model.',
  'mqmcofactors' =>
    'List of cofactors to be analysed in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}}.',
  'phenocol' =>
    'Column number in the phenotype matrix which should be used as the phenotype. This can be a vector of integers.',
  'verbose' =>
    'Display more output on verbose=TRUE'
}

ARGV.each do | fn | 

  raise "File not found #{fn}!" if !File.exist?(fn)
  print "\nParsing #{fn}..."

  buf = nil
  File.open(fn) { | f | buf = f.read }
  
  # parse buffer and strip between inputs
  outbuf = []
  inside_input = false
  skipone = false
  buf.each do | s |
    if s.strip =~ /^%\s+\\input\{\"(\S+?)\"\}/
      inputfn = $1
      print "\nInjecting #{inputfn}"
      inside_input = true
      # inject inputfn
      raise "File not found #{inputfn}!" if !File.exist?(inputfn)
      outbuf.push s
      outbuf.push File.new(inputfn).read
      outbuf.push "% -----^^ "+inputfn.strip+" ^^-----\n"
    end
    # Now check keywords
    REPL.each do | k, v |
      if s.strip =~ /%\s+\\#{k}/
        outbuf.push v + " % \\"+k+"\n"
        skipone = true
        next
      end
    end
    inside_input = false if s.strip == '}'
    outbuf.push s if !inside_input and !skipone
    skipone = false
  end
  File.open(fn,"w") do | f |
    f.print outbuf
  end

end
