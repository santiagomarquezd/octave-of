git status | grep modified | cut -d: -f 2 | xargs git add
