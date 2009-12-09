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
# closing curly brace.
#
# Usage:
#
#   ruby ./contrib/scripts/repl_inputs.rb inputfile(s)
#
# Example
#
#   ruby ./contrib/scripts/repl_inputs.rb man/*.Rd
#

ARGV.each do | fn | 

  raise "File not found #{fn}!" if !File.exist?(fn)
  print "\nParsing #{fn}..."

  buf = nil
  File.open(fn) { | f | buf = f.read }
  
  # parse buffer and strip between inputs
  outbuf = []
  input = false
  buf.each do | s |
    if s.strip =~ /^%\s+\\input\{\"(\S+?)\"\}/
      inputfn = $1
      input = true
      # inject inputfn
      raise "File not found #{inputfn}!" if !File.exist?(inputfn)
      outbuf.push s
      outbuf.push File.new(inputfn).read
      outbuf.push "% -----^^ "+inputfn.strip+" ^^-----\n"
    end
    input = false if input and s.strip == '}'
    outbuf.push s if !input
  end
  File.open(fn,"w") do | f |
    f.print outbuf
  end

end
