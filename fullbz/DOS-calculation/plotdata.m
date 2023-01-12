subcalcname = "sg81-9-3098-4990";
subcalcname = "sg81-10-4822-4358";
subcalcname = "sg13-9-1916-4417";
subcalcname = "sg81-9-3098-4384";
dosdata = importdata(subcalcname+ "-output.txt");
figure;
plot(dosdata(2:end, 1), dosdata(2:end, 2))
% ylim([0, 80]);
xlabel("Frequency");
ylabel("DOS");
