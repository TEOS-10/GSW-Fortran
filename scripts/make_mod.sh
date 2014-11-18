#!/bin/sh
#
# Create module interface statements from GSW .f90 files.
#
# Usage: ./make_mod.sh *.f90
#

toolbox_name="gsw_mod_toolbox"

outfile=${toolbox_name}".f90"

tmpfile="__temp.$$"

echo "module "${toolbox_name} > $outfile
cat << END >> $outfile

implicit none

END

while [ -n "$1" ]; do

	mod=`grep "module " $1`
	if [ -n "$mod" ]; then
		echo "ignoring module: "$1
		shift
		continue
	fi
	mod=`grep "elemental " $1`
	if [ -z "$mod" ]; then
		mod=`grep "pure " $1`
		if [ -z "$mod" ]; then
			echo "ignoring non-elemental/pure routine: "$1
			shift
			continue
		fi
	fi

	echo "public :: "`basename $1 .f90` >> $outfile

	sed -n '
/^elemental sub.*& *$/ {
	N
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^elemental func.*& *$/ {
	N
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^pure sub.*& *$/ {
	N
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^pure func.*& *$/ {
	N
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^elemental func/ {
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^elemental sub/ {
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^pure func/ {
	p
	s/^.*\(gsw_[a-z0-9_]*\) *(.*/\1/
	h
	b
}
/^pure sub/ {
	p
	b
}
/:: *gsw_.*/ {
	p
	b
}
/^end func/ {
	G
	s/\n/ /
	p
	b
}
/^end sub/ {
	G
	s/\n/ /
	p
	b
}
/selected_real/ { p }
/implicit / { p }
/intent(in)/ { p }
/intent(out)/ { p }
' $1 >> $tmpfile

	echo "" >> $tmpfile
	shift
done

cat << END >> $outfile

interface

END

sed 's/\(.*\)/    \1/' $tmpfile >> $outfile
/bin/rm $tmpfile

cat << END >> $outfile
end interface

END
echo "end module "${toolbox_name} >> $outfile
