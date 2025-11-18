#!/bin/bash
make
./wsd example/test_ref.fa -m example/test_tems.fa > example/test_user_test_out.txt

if cmp --silent "example/test_user_test_out.txt" "example/test_ref_out.txt"; then
    echo "WSD test successfully!"
else
    echo "WSD test failed!"
fi