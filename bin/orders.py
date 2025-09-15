str1 = "/home/uab/uab513084/trait-selection/test-2-directional/bin/run_all_AL.sh "
file = open("orders.dat","a")#write mode
for i in range(1,112):
    a = str(i).zfill(5)
    str2 = "/home/uab/uab513084/trait-selection/test-2-directional/bin/running/{number}/".format(number = a)
    str3 = " /home/uab/uab513084/trait-selection/test-2-directional/bin 10 {number}".format(number = a)
    string = str1 + str2 + str3 + "\n"
    file.write(string)
str4 = "/home/uab/uab513084/trait-selection/test-2-directional/bin/subs.sh /home/uab/uab513084/trait-selection/test-2-directional/bin 300000 1000"
file.write(str4)
file.close()
