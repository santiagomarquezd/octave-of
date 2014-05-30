git status | grep deleted | cut -d ":" -f 2 | xargs git checkout --

