
#!/bin/bash

cd ../build
make
if [ $? -ne 0 ] ; then
    echo "";
    echo "error in the build" ;
    echo "";
else 
    cd ../run
    ./odt.x
fi

