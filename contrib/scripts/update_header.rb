#! /usr/bin/ruby
#
# This script updates the header of a source file
#
# by Pjotr Prins

require 'parsedate'

def git_info(source, mainauthor)
  git = `git log #{source}`.split(/\n/)
  git[2] =~ /Date: /
  lastmodified = ParseDate.parsedate($'.strip)
  modifiedby = []
  git.grep(/^Author/).uniq.each do | author |
    if author !~ /#{mainauthor}/
      author =~ /Author: (.*) </
      modifiedby.push $1
    end
  end
  return lastmodified, modifiedby
end

def file_R(buf, source)
  outbuf = []
  # parse for modified by
  lastmodified, modifiedby = git_info(source,"anny")
  d = lastmodified
  t = Time.mktime(d[0],d[1],d[2])
  # parse for methods
  methods = ""
  buf.each do | s |
    if s =~ /^(\S+)\s+<-/
      methods += $1 + "\n#           "
    end
  end
  inside_header = true
  buf.each do | s |
    if inside_header and s !~ /^#/
      # now inject header
      outbuf.push <<R_HEADER
#####################################################################
#
# #{File.basename(source)}
#
# Copyright (c) 2009, Danny Arends
#
# Modified by #{modifiedby.join(" and ")}
#
# 
# first written Februari 2009
# last modified #{t.strftime("%B %Y")}
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
# Part of the R/qtl package
# Contains: #{methods}
#
#####################################################################


R_HEADER
      inside_header = false
    else
      outbuf.push s if !inside_header
    end
  end
  outbuf
end

def file_C(buf, source)
  outbuf = []
  # parse for modified by
  lastmodified, modifiedby = git_info(source,"ansen")
  d = lastmodified
  t = Time.mktime(d[0],d[1],d[2])
  inside_header = true
  buf.each do | s |
    if inside_header and s !~ /^[\s\/]*\*/
      # now inject header
      outbuf.push <<C_HEADER
/**********************************************************************
 *
 * #{File.basename(source)}
 *
 * Copyright (c) 2009 Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * Modified by #{modifiedby.join(" and ")}
 *
 * first written before 2000
 * last modified #{t.strftime("%B %Y")}
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 *
 **********************************************************************/

C_HEADER
      inside_header = false
    else
      outbuf.push s if !inside_header
    end
  end
  outbuf
end


ARGV.each do | fn | 

  raise "File not found #{fn}!" if !File.exist?(fn)
  print "\nParsing #{fn}..."

  buf = nil
  File.open(fn) { | f | buf = f.read }
  
  # parse buffer and strip header replacing it with new

  if fn =~ /\.R/
    outbuf = file_R(buf, fn)
  elsif fn =~ /\.[cChH]/
    outbuf = file_C(buf, fn)
  else
    raise 'Unknown file extension for '+fn
  end

  File.open(fn,"w") do | f |
    f.print outbuf
  end

end
