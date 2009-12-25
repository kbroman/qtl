#! /usr/bin/ruby
#
# This script updates the header of a source file
#
# by Pjotr Prins



R_HEADER =<<R_BLOCK
#####################################################################
#
# mqmaugment.R
#
# copyright (c) 2009, Danny Arends, Pjotr Prins and Karl W. Broman
# last modified July, 2009
# first written Feb, 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
#####################################################################
R_BLOCK


def file_R(buf)
  outbuf = []
  inside_header = false
  buf.each do | s |
    if inside_header and s !~ /^#/
      # now inject header
      outbuf.push "# new header!\n"
    else
      outbuf.push s
    end
  end
end


ARGV.each do | fn | 

  raise "File not found #{fn}!" if !File.exist?(fn)
  print "\nParsing #{fn}..."

  buf = nil
  File.open(fn) { | f | buf = f.read }
  
  # parse buffer and strip header replacing it with new

  if fn =~ /\.R/
    outbuf = file_R(buf)
  else
    raise 'Unknown file extension for '+fn
  end

  File.open(fn,"w") do | f |
    f.print outbuf
  end

end
