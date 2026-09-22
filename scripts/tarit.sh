#!/bin/zsh
# create a tar archive

# deduce version automatically from the appropriate include file
version=`grep 'Version' README | head -1 | sed 's/.*Version *//' | sed "s/ *\*.*//"`

#echo $version
#exit


origdir=`pwd | sed 's/.*\///'`
echo "Will make an archive of $origdir/"

dirhere=`pwd`
dirhere=$dirhere:t

dirtar=JetVHeto-$version
tarname=$dirtar.tgz
tmptarname=tmp-$tarname

# make sure we have Makefile with use CGAL=no

if [[ -e ../$tarname ]]
then
  echo "Tarfile $tarname already exists. Not proceeding."
elif [[ -e /tmp/$dirtar ]]
then
  echo "/tmp/$dirtar already exists, not proceeding"
else
  pushd ..

  echo "Creating tmp-$tarname"
  tar -h --exclude '.svn*' --exclude '*~' \
      -zcf $tmptarname \
                      $dirhere/(src|scripts)/*.(f90|sh|pl) \
                      $dirhere/ChangeLog \
                      $dirhere/Makefile\
                      $dirhere/(AUTHORS|BUGS|COPYING|README|NEWS)\
                      $dirhere/modules/.dummy\
                      $dirhere/obj/.dummy\
		      $dirhere/fixed-order/*.fxd \
		      $dirhere/fixed-order-masses/*.fxd \
                      $dirhere/n3lo-fixed-order/*.fxd \
                      $dirhere/benchmarks-n3lo-smallR/*.res \
		      $dirhere/python/*   

  
  fulltarloc=`pwd`
  pushd /tmp
  echo "Unpacking it as /tmp/$dirhere"
  tar zxf $fulltarloc/$tmptarname
  mv -v /tmp/$dirhere /tmp/$dirtar
  echo "Repacking it with directory name $dirtar"
  tar zcvf $fulltarloc/$tarname $dirtar 
  echo 
  echo "Removing /tmp/$dirhere"
  rm -rf $dirtar
  popd
  rm -v $tmptarname

  echo ""
    # if it's gavin running this then automatically copy the tarfile
    # to the web-space
  # webdir=~salam/www/repository/software/fastjet/
  # if [[ $USER = salam && -e $webdir ]]
  #     then
  #     echo "Copying .tgz file to web-site"
  #     cp -vp $tarname $webdir
  #     echo "************   Remember to edit web page **********"
  # fi

  popd

  # reminders about what to do for svn
  URL=`svn info | grep URL | sed 's/^.*URL: //'`
  echo $URL
  tagURL=`echo $URL | sed "s/trunk.*JetVHeto/tags\/JetVHeto-$version/"`
  echo "Remember to tag the version:"
  echo "   svn copy  -m 'tagged release of release $version' $URL $tagURL"
  echo "Copy it to hepforge:"
  echo "   scp -p $fulltarloc/$tarname login.hepforge.org:/hepforge/home/jetvheto/downloads/"
  echo "Then edit web page on hepforge"
fi

#tar zcf $tarname


