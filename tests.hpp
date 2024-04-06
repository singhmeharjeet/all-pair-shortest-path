#include <stdio.h>

void check(bool res, char* name) {
	if (res) {
		printf("Passed");
	} else {
		printf(">>----> Failed");
	}
	printf(" %s\n", name);
}

bool test_search_pass() {
	return false;
}
