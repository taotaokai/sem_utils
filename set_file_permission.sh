#!/bin/bash

find utils -type f | xargs chmod 644
find utils -type f -name "*.py" | xargs chmod 755 
find utils -type f -name "*.sh" | xargs chmod 755 