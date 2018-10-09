exp = readExp();
exp = genExp(exp);

subjid = "1" + Math.random().toString().substring(3,8);
file_name = subj_id + ".csv";

stage = "train";
trial_idx = -1;
in_trial = false;

RTs = [];
path = [];
cur = -1;
start = -1;
goal = -1;

nextTrial();


// Read experiment template from textarea
//
function readExp() {
    var exp = {};
    var lines = $("#experiment").val().split("\n");
    exp.N = parseInt(lines[0], 10);

    // read adjacency 
    exp.adj = [];
    for (var i = 1; i <= exp.N; i++) {
        var a = lines[i].trim().split(" ");
        var b = [];
        for (var j = 0; j < 4; j++) {
            b.push(parseInt(a[j], 10));
        }
        exp.adj.push(b);
    }
    l = exp.N + 1; // current line

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

    // read state names
    exp.names = [];
    for (var i = 0; i < exp.N; i++) {
        exp.names.push(lines[l].trim());
        l++;
    }

    return exp;
}


function genExp(exp) {
    // shuffle state names
    exp.names.sort(function(a, b) {return 0.5 - Math.random()});

    // generate training trials
    exp.train_trials = genTrials(exp.train);

    // generate test trials
    exp.test_trials = genTrials(exp.test);

    // optionally rotate graph TODO enable
    /*
    exp.rotate = Math.floor(Math.random() * 4);
    for (var i = 0; i <= exp.N; i++) {
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


function nextTrial() {
    $("#trial_page").hide();
    $("#message").text("");
    in_trial = false;

    trial_idx++;
    if (stage == "train") {
        trials = exp.train_trials;
    } else {
        trials = exp.test_trials;
    }

    if trial_idx >= trials.length) {
        if (stage == "train") {
            // kick off test phase
            stage = "test";
            trial_idx = -1;
            $("#test_page").show();
            sleep(5000).then(() => {
                $("#test_page").hide();
                nextTrial();
            });
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
        sleep(100).then(() => {
            $("#countdown").text("2...");
            sleep(100).then(() => {
                $("#countdown").text("1...");
                sleep(100).then(() => {
                    $("#countdown").text("GO!");
                    in_trial = true;
                    last_keypress_time = Date.getTime();
                    sleep(100).then(() => {
                        $("#countdown").text("");
                    });
                });
            });
        });
    });
}


function checkKeyPressed(e) {
    e = e || window.event;
    console.log("key press" + e.which);

    if (in_trial) {
        RT = Date.getTime() - last_keypress_time;
        last_keypress_time = Date.getTime();
        RTs.push(RT);
        RT_tot += RT;
        var next = cur;
        $("#message").text("");

        if ((e).keyCode == "39") {
            next = exp.adj[cur][0];
        } else if ((e).keyCode == "38") {
            next = exp.adj[cur][1];
        } else if ((e).keyCode == "37") {
            next = exp.adj[cur][2];
        } else if ((e).keyCode == "40") {
            next = exp.adj[cur][3];
        }

        if (next >= 0) {
            cur = next;
            redraw();
        }

        if ((e).key === ' ' || (e).key === 'Spacebar') {
            if (cur == goal) {
                $("#message").css("color", "green");
                $("#message").text("SUCCESS!!");
                logTrial();
                sleep(1000).then(() => {
                    nextTrial();
                });
            } else {
                $("#message").css("color", "red");
                $("#message").text("Incorrect");
            }
        }
    }

    // at some poin t
    // in_trial = false;
}


function logTrial() {
    RT_str = (RTs.toString()).replace(/,/g/," ");
    path_str = (path.toString()).replace(/,/g/," ");
    console.log(path_str);
    console.log(RT_str);
    row = "A," + subj_id + "," + stage + "," + start.toString() + "," + goal.toString() + "," + path_str + "," + RT_str + "," + RT_tot.toString() + "\n";
    $.post("results_data.php", {postresult: row, postfile: file_name});
}


function redraw() {
    cur_name = exp.names[cur - 1];
    start_name = exp.names[start - 1];
    goal_name = exp.names[goal - 1];
    var adj_names = [];
    for (var i = 0; i < 4; i++) {
        if (adj[cur][i] <= 0) {
            adj_names.push("&#11044;");
        } else {
            adj_names.push(exp.names[exp.adj[cur][j]] - 1);
        }
    }

    $("#cur_state").text(cur_name);
    $("#goal_state").text("Go to " + goal_name);
    $("#right_state").text(adj_names[0]);
    $("#up_state").text(adj_names[1]);
    $("#left_state").text(adj_names[2]);
    $("#down_state").text(adj_names[3]);

    $("#from_state").text(start_name);
    $("#to_state").text(goal_name);
}


// helper f'n
function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

