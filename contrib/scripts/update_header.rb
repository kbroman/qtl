#! /usr/bin/ruby
#
# This script updates the header of a source file
#
# by Pjotr Prins

require 'parsedate'

class GitLogEntry
  attr_reader :id, :author, :date, :comment
  def initialize buf
    raise "Git log problem "+buf.to_s if buf[0] !~ /^commit/
    @id = buf[0]
    pos = 1
    pos += 1 if buf[pos] =~ /^Merge:/
    buf[pos] =~ /^Author: (.*) </
    @author = $1.strip
    @author = 'Danny Arends' if @author == 'DannyArends'
    pos += 1
    buf[pos] =~ /^Date:   /
    @date = $'.strip
    @comment = buf[3..-1]
  end
end

class GitLog
  def initialize fn
    @list = []
    buf = []
    # print "Fetching git log from #{fn}\n"
    `git log #{fn}`.split(/\n/).each do | s |
      if s =~ /^commit/ and buf.size>0
        @list.push GitLogEntry.new(buf)
        buf = []
      end
      buf.push s
    end
    @list.push GitLogEntry.new(buf)
  end

  def modified
    @list.each do | commit |
      next if commit.comment.join =~ /header/i
      return ParseDate.parsedate(commit.date)
    end
    'unknown'
  end

  def authors
    modifiedby = []
    @list.each do | commit |
      next if commit.comment.join =~ /header/i
      modifiedby.push commit.author
    end
    modifiedby.uniq
  end
end

def git_info(source, mainauthor)
  gitlog = GitLog.new(source)
  modifiedby = []
  return gitlog.modified, gitlog.authors
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
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by #{modifiedby.join(" and ")}
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
  elsif fn =~ /\.o$/
    next
  else
    raise 'Unknown file extension for '+fn
  end

  File.open(fn,"w") do | f |
    f.print outbuf
  end

end
