function readTextFile(file)
{
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", file, false);
    rawFile.onreadystatechange = function ()
    {
        if(rawFile.readyState === 4)
        {
            if(rawFile.status === 200 || rawFile.status == 0)
            {
                var allText = rawFile.responseText;
                alert(allText);
            }
        }
    }
    rawFile.send(null);
}

function foo()
{
    document.getElementById("apple").innerHTML = "Hey, Bro!";
}

function readFile()
{
    if (window.File && window.FileReader && window.FileList && window.Blob) {
        var fname = "file:///Users/moorer/projects/ZetaHunter/output/degapped_alignment.fa"
        var reader = new FileReader();

        reader.onload = function(e) {
            document.getElementById("pie").innerHTML = reader.result;
        }

        reader.readAsText(fname);

    } else {
          alert('The File APIs are not fully supported by your browser.');
    }
}
