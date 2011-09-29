# check
cd ~/R
R --vanilla CMD check qcc

# build
cd ~/R
R --vanilla CMD build --force qcc

# install
lib="/Library/Frameworks/R.framework/Resources/library"
sudo R --vanilla CMD REMOVE --library="$lib" qcc
sudo R --vanilla CMD INSTALL --library="$lib" qcc_2.2.tar.gz # <-- change with the last version

# create by hand a zip archive for MS Windows
# lib="/Library/Frameworks/R.framework/Resources/library"
# cd ~/tmp
# cp $lib/qcc .
# zip -r qcc_2.2.zip qcc
# mv qcc_2.2.zip ~/R
# rm -rf qcc
