#!/bin/bash

rm -rf package
rm package_macOS_10_15_5_x86_64bit.zip
cp -av ./package_raw ./package

pyinstaller substructure_search.py
mv ./dist/substructure_search ./package/service
mkdir ./package/service/pygly/
cp ~/codes/python2/PyGly/pygly/*.ini ./package/service/pygly/
# TODO maybe also include glytoucan version?
cp ./glycan_set_glygen.tsv ./package/service/glycan_set_glygen.tsv
cp ./index.html ./package/service/index.html
cp ./service.ini ./package/service/service.ini
rm -rf build
rm -rf dist
rm substructure_search.spec

pyinstaller submit.py
mv ./dist/submit ./package/submit
chmod +x ./package/submit/submit
rm -rf build
rm -rf dist
rm submit.spec

chmod +x ./package/service/substructure_search


zip -r package_macOS_10_15_5_x86_64bit.zip package

