document.addEventListener("DOMContentLoaded", init, false);


function init() {
    var aliveCells = [];
    var gridSize = 5;
    var gridHtml = createGrid(gridSize);
    document.querySelector('#grid').innerHTML = gridHtml;
    mapEventToCells();
    aliveCells = [0, 11, 23, 26, 25, 29, 30, 31];
    displayAlive(aliveCells);

    function mapEventToCells() {
        var cells = document.querySelectorAll(".cell");
        for (let i = 0; i < cells.length; i++) {
            const element = cells[i];
            element.addEventListener("click", function () {
                this.classList.toggle("alive");
                if(this.classList.value.indexOf("alive") !== -1) {
                    aliveCells.push(1)
                }
                
            })
        }
    }

    function createGrid(size) {
        let gridHtml = '';
        let numberOfCells = size*size;
        for (let i = 0; i < numberOfCells; i++) {
            gridHtml += '<dir class="cell"></dir> \n';
        }
        return gridHtml;
    }
    function displayAlive(aliveCells) {
        var cells = document.querySelectorAll(".cell");
        for (let i = 0; i < aliveCells.length; i++) {
            cells[aliveCells[i]].classList.add("alive");
        }
    }

}

function testFunctions(){
    
}