in_trial = false;

function initExp() {
    console.log("initExp");

    exp = readExp();
    exp = genExp(exp);

    subj_id = "1" + Math.random().toString().substring(3,8);
    file_name = 'results/' + subj_id + ".csv";

    stage = "train";
    trial_idx = -1;
    in_trial = false;

    RTs = [];
    keys = [];
    path = [];
    cur = -1;
    start = -1;
    goal = -1;

    $.post("results_data.php", {postresult: "group, subj_id, stage, start, goal, path, length, RTs, keys, RT_tot, timestamp, datetime\n", postfile: file_name })

    nextTrial();
}


// Read experiment template from textarea
//
function readExp() {
    console.log("readExp");

    var exp = {};
    var lines = $("#experiment").val().split("\n");
    exp.N = parseInt(lines[0], 10);

    // read adjacency 
    exp.adj = [];
    exp.x = [];
    exp.y = [];
    for (var i = 1; i <= exp.N; i++) {
        var a = lines[i].trim().split(" ");
        var b = [];
        for (var j = 0; j < 4; j++) {
            b.push(parseInt(a[j], 10));
        }
        exp.x.push(parseInt(a[4], 10));
        exp.y.push(parseInt(a[5], 10));
        exp.adj.push(b);
    }
    l = exp.N + 1; // current line

    // read state names
    exp.names = [];
    for (var i = 0; i < exp.N; i++) {
        exp.names.push(lines[l].trim());
        l++;
    }

    // read training tasks
    exp.ntrain = parseInt(lines[l], 10);
    l++;
    exp.train = [];
    for (var i = 0; i < exp.ntrain; i++) {
        var a = lines[l].trim().split(" ");
        l++;
        var task = {};
        task.s = parseInt(a[0], 10);
        task.g = parseInt(a[1], 10);
        task.n = parseInt(a[2], 10);
        exp.train.push(task);
    }

    // read test tasks
    exp.ntest = parseInt(lines[l], 10);
    l++;
    exp.test = [];
    for (var i = 0; i < exp.ntest; i++) {
        var a = lines[l].trim().split(" ");
        l++;
        var task = {};
        task.s = parseInt(a[0], 10);
        task.g = parseInt(a[1], 10);
        task.n = parseInt(a[2], 10);
        exp.test.push(task);
    }

    return exp;
}


function genExp(exp) {
    console.log("genExp");

    // shuffle state names
    exp.names.sort(function(a, b) {return 0.5 - Math.random()});

    // generate training trials
    exp.train_trials = genTrials(exp.train);

    // generate test trials
    exp.test_trials = genTrials(exp.test);

    // optionally rotate graph
    /*
    DON'T DO! it doesn't rotate the full map; just messes up stuff
    exp.rotate = Math.floor(Math.random() * 4);
    for (var i = 0; i < exp.N; i++) {
        var a = exp.adj[i].slice();
        for (var j = 0; j < 4; j++) {
            exp.adj[i][j] = a[(j + exp.rotate) % 4];
        }
    }
    */
    return exp;
}


function genTrials(desc) {
    trials = [];
    for (var i = 0; i < desc.length; i++) {
        for (var j = 0; j < desc[i].n; j++) {
            var task = {};
            task.s = desc[i].s;
            task.g = desc[i].g;
            while (task.s <= 0) {
                task.s = Math.floor(Math.random() * exp.N) + 1;
            }
            while (task.g <= 0 || task.g == task.s) {
                task.g = Math.floor(Math.random() * exp.N) + 1;
            }
            trials.push(task);
        }
    }
    trials.sort(function(a, b) {return 0.5 - Math.random()});
    return trials;
}

// Fisher-Yates (aka Knuth) Shuffle, from 
// https://stackoverflow.com/questions/2450954/how-to-randomize-shuffle-a-javascript-array
//
function shuffle(array) {
  var currentIndex = array.length, temporaryValue, randomIndex;

  // While there remain elements to shuffle...
  while (0 !== currentIndex) {

    // Pick a remaining element...
    randomIndex = Math.floor(Math.random() * currentIndex);
    currentIndex -= 1;

    // And swap it with the current element.
    temporaryValue = array[currentIndex];
    array[currentIndex] = array[randomIndex];
    array[randomIndex] = temporaryValue;
  }

  return array;
}


function nextTrial() {
    console.log("nextTrial " + trial_idx);

    $("#trial_page").hide();
    $("#message").text("");
    in_trial = false;

    trial_idx++;
    if (stage == "train") {
        trials = exp.train_trials;
    } else {
        trials = exp.test_trials;
    }

    if (trial_idx >= trials.length) {
        if (stage == "train") {
            // kick off test phase
            stage = "test";
            trial_idx = -1;
            $("#test_page").show();
        } else {
            // fin
            $("#final_page").show();
        }
        return;
    }

    start = trials[trial_idx].s;
    cur = start;
    goal = trials[trial_idx].g;

    RT_tot = 0;
    RTs = [];
    keys = [];
    path = [cur];

    redraw();
    $("#new_trial_page").show();

    // countdown
    sleep(2000).then(() => {
        $("#new_trial_page").hide();
        $("#trial_page").show();
        $("#countdown").text("3...");
        stateColor("grey");
        sleep(1000).then(() => {
            $("#countdown").text("2...");
            sleep(1000).then(() => {
                $("#countdown").text("1...");
                sleep(1000).then(() => {
                    $("#countdown").text("GO!");
                    stateColor("white");
                    in_trial = true;
                    redraw();
                    last_keypress_time = (new Date()).getTime();
                    sleep(1000).then(() => {
                        $("#countdown").text("");
                    });
                });
            });
        });
    });
}


function stateColor(color) {
    $("#cur_state").css("color", color);
    $("#right_state").css("color", color);
    $("#up_state").css("color", color);
    $("#left_state").css("color", color);
    $("#down_state").css("color", color);
}


function checkKeyPressed(e) {
    var e = window.event || e;

    if (in_trial) {
        console.log("key press " + e.which);

        RT = (new Date()).getTime() - last_keypress_time;
        last_keypress_time = (new Date()).getTime();
        RTs.push(RT);
        keys.push((e).keyCode);
        RT_tot += RT;
        var next = -1;
        $("#message").text("");

        // get next state
        if ((e).keyCode == "39") {
            next = exp.adj[cur - 1][0];
        } else if ((e).keyCode == "38") {
            next = exp.adj[cur - 1][1];
        } else if ((e).keyCode == "37") {
            next = exp.adj[cur - 1][2];
        } else if ((e).keyCode == "40") {
            next = exp.adj[cur - 1][3];
        }

        if (stage == "train") {
            // move to next state 
            if (next >= 0) {
                cur = next;
                stateColor("grey");
                in_trial = false;
                path.push(next);
               // sleep(750).then(() => {
                    stateColor("white");
                    in_trial = true;
                    redraw();
               // });
            }

            // if goal is reached => start next trial
            if ((e).key === ' ' || (e).key === 'Spacebar') {
                if (cur == goal) {
                    $("#message").css("color", "green");
                    $("#message").text("SUCCESS!!");
                    in_trial = false;
                    logTrial();
                    sleep(1000).then(() => {
                        nextTrial();
                    });
                } else {
                    $("#message").css("color", "red");
                    $("#message").text("Incorrect");
                }
            }
        } else { // stage == "test"
            // end trial after first button press
            if (next >= 0) {
                path.push(next);
                stateColor("grey");
                in_trial = false;
                redraw();
                logTrial();
                sleep(1000).then(() => {
                    stateColor("white");
                    nextTrial();
                });
            }
        }
    }

    return true;
}


function logTrial() {
    var RT_str = (RTs.toString()).replace(/,/g, ' ');
    var path_str = (path.toString()).replace(/,/g, ' ');
    var key_str = (keys.toString()).replace(/,/g, ' ');
    var d = new Date();
    var t = d.getTime() / 1000;
    var row = "A," + subj_id + "," + stage + "," + start.toString() + "," + goal.toString() + "," + path_str + "," + path.length.toString() + "," + RT_str + "," + key_str + "," + RT_tot.toString() + "," + t.toString() + "," + d.toString() + "\n";
    console.log(row);
    $.post("results_data.php", {postresult: row, postfile: file_name});
}


function redraw() { 
    cur_name = exp.names[cur - 1];
    start_name = exp.names[start - 1];
    goal_name = exp.names[goal - 1];
    if (in_trial) {
        white = "white";
        green = "green";
    } else {
        white = "grey";
        green = "#008800";
    }

    var canvas = document.getElementById("map");
    var ctx = canvas.getContext("2d");
    ctx.fillStyle = "black";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    scale = 300;

    xoffs = canvas.width / 2 - (Math.max(...exp.x) - Math.min(...exp.x)) * scale / 2;
    yoffs = canvas.height / 2 - (Math.max(...exp.y) - Math.min(...exp.y)) * scale / 2;

    // edges
    for (var i = 0; i < exp.N; i++) {
        for (var j = 0; j < 4; j++) {
            next = exp.adj[i][j] - 1;
            if (next > 0) {
                ctx.beginPath();
                ctx.strokeStyle = white;
                ctx.lineWidth = 5;
                ctx.moveTo(exp.x[i] * scale + xoffs, exp.y[i] * scale + yoffs);
                ctx.lineTo(exp.x[next] * scale + xoffs, exp.y[next] * scale + yoffs);
                ctx.stroke();
            }
        }
    }

    // nodes
    for (var i = 0; i < exp.N; i++) {
        ctx.beginPath();
        ctx.rect(exp.x[i] * scale + xoffs - 0.25 * scale, exp.y[i] * scale + yoffs - 0.25 * scale, 0.5 * scale, 0.50 * scale);
        ctx.fillStyle = "black"; 
        ctx.fill();
        if (i == cur - 1) {
            ctx.lineWidth = 15;
        } else {
            ctx.lineWidth = 2;
        }
        ctx.strokeStyle = white;
        ctx.stroke();
        ctx.lineWidth = 0;

        if (i == goal - 1) {
            ctx.fillStyle = green;
            ctx.font = "bold 25px Ariel";
        } else {
            ctx.fillStyle = white;
            ctx.font = "25px Ariel";
        }
        ctx.textAlign = "center";
        ctx.fillText(exp.names[i], exp.x[i] * scale + xoffs, exp.y[i] * scale + yoffs);
    }

    $("#goal_state").text("Go to " + goal_name);

    $("#from_state").text(start_name);
    $("#to_state").text(goal_name);
}


// helper f'n
function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

