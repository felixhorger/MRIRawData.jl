MRIRawData.jl
-------------

Reading common MRI manufacturer's raw data formats.

Currently supported is only Siemens twix.
Please do contribute!

Requirements
------------
```
cd whereeveryouwanttohaveit
git clone https://github.com/felixhorger/pymapvbvd # My fork removing the printing of progress
cd pymapvbvd
python3 setup.py install --user
```

Credits
-------
For reading Siemens data, this wraps [pymapvbvd](https://github.com/wtclarke/pymapvbvd) by William Clarke.
I can't tell you how grateful I am somebody took the time to port it, making the MATLAB version obsolete,
cheers to that!

